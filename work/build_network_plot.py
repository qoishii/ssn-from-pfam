#!/usr/bin/env python3
"""Build and plot a sequence similarity network from BLAST-like hits."""

from __future__ import annotations

import argparse
import csv
import math
import random
import struct
import subprocess
import sys
import xml.etree.ElementTree as ET
import zlib
from collections import defaultdict, deque
from pathlib import Path


try:
    import networkx as nx
except ImportError:
    nx = None


PALETTE = [
    (31, 119, 180),
    (255, 127, 14),
    (44, 160, 44),
    (214, 39, 40),
    (148, 103, 189),
    (140, 86, 75),
    (227, 119, 194),
    (127, 127, 127),
    (188, 189, 34),
    (23, 190, 207),
]

HEX_PALETTE = [f"#{red:02x}{green:02x}{blue:02x}" for red, green, blue in PALETTE]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build an SSN edge table and render an SVG network figure."
    )
    parser.add_argument("--similarity", type=Path, required=True, help="Input similarity.tsv.")
    parser.add_argument("--metadata", type=Path, required=True, help="Filtered metadata TSV.")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory.")
    parser.add_argument("--min-identity", type=float, default=50.0, help="Minimum percent identity.")
    parser.add_argument("--min-coverage", type=float, default=0.70, help="Minimum min(query, target) coverage.")
    parser.add_argument("--max-evalue", type=float, default=1e-5, help="Maximum e-value.")
    parser.add_argument("--width", type=int, default=1800, help="PNG width in pixels.")
    parser.add_argument("--height", type=int, default=1400, help="PNG height in pixels.")
    parser.add_argument("--seed", type=int, default=7, help="Deterministic layout seed.")
    parser.add_argument("--label-font-size", type=int, default=18, help="Node label font size for Graphviz renders.")
    parser.add_argument(
        "--legend-components",
        type=int,
        default=8,
        help="Number of largest components to list in the Graphviz legend.",
    )
    return parser.parse_args()


def read_metadata(path: Path) -> dict[str, dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["Entry"]: row for row in reader}


def accession_from_blast_id(identifier: str) -> str:
    parts = identifier.split("|")
    if len(parts) >= 2:
        return parts[1]
    return identifier


def read_edges(
    path: Path, min_identity: float, min_coverage: float, max_evalue: float
) -> dict[tuple[str, str], dict[str, float | str]]:
    best_edges: dict[tuple[str, str], dict[str, float | str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            query = accession_from_blast_id(row["query"])
            target = accession_from_blast_id(row["target"])
            if query == target:
                continue
            pident = float(row["pident"])
            coverage = float(row["coverage"])
            evalue = float(row["evalue"])
            bitscore = float(row["bitscore"])
            score = float(row["score"])
            if pident < min_identity or coverage < min_coverage or evalue > max_evalue:
                continue
            source, dest = sorted([query, target])
            key = (source, dest)
            existing = best_edges.get(key)
            if existing is None or bitscore > float(existing["bitscore"]):
                best_edges[key] = {
                    "source": source,
                    "target": dest,
                    "pident": pident,
                    "coverage": coverage,
                    "evalue": evalue,
                    "score": score,
                    "bitscore": bitscore,
                }
    return best_edges


def connected_components(nodes: list[str], edges: dict[tuple[str, str], dict[str, float | str]]) -> dict[str, int]:
    adjacency: dict[str, list[str]] = {node: [] for node in nodes}
    for source, target in edges:
        adjacency[source].append(target)
        adjacency[target].append(source)

    components: list[list[str]] = []
    seen: set[str] = set()
    for node in nodes:
        if node in seen:
            continue
        queue = deque([node])
        seen.add(node)
        component = []
        while queue:
            current = queue.popleft()
            component.append(current)
            for neighbor in adjacency[current]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    queue.append(neighbor)
        components.append(component)

    components.sort(key=lambda item: (-len(item), item[0]))
    cluster_by_node: dict[str, int] = {}
    for index, component in enumerate(components, start=1):
        for node in component:
            cluster_by_node[node] = index
    return cluster_by_node


def degree_by_node(nodes: list[str], edges: dict[tuple[str, str], dict[str, float | str]]) -> dict[str, int]:
    degrees = {node: 0 for node in nodes}
    for source, target in edges:
        degrees[source] += 1
        degrees[target] += 1
    return degrees


def build_networkx_graph(
    nodes: list[str],
    metadata: dict[str, dict[str, str]],
    edges: dict[tuple[str, str], dict[str, float | str]],
):
    if nx is None:
        return None
    graph = nx.Graph()
    for node in nodes:
        row = metadata[node]
        graph.add_node(
            node,
            protein_name=row.get("Protein names", ""),
            organism=row.get("Organism", ""),
            length=row.get("Length", ""),
            reviewed=row.get("Reviewed", ""),
            ec_number=row.get("EC number", ""),
        )
    for edge in edges.values():
        graph.add_edge(
            str(edge["source"]),
            str(edge["target"]),
            pident=float(edge["pident"]),
            coverage=float(edge["coverage"]),
            evalue=float(edge["evalue"]),
            score=float(edge["score"]),
            bitscore=float(edge["bitscore"]),
        )
    return graph


def networkx_components_and_degrees(graph) -> tuple[dict[str, int], dict[str, int]]:
    components = sorted(nx.connected_components(graph), key=lambda item: (-len(item), sorted(item)[0]))
    clusters: dict[str, int] = {}
    for index, component in enumerate(components, start=1):
        for node in component:
            clusters[node] = index
    degrees = {node: int(degree) for node, degree in graph.degree()}
    nx.set_node_attributes(graph, clusters, "cluster")
    nx.set_node_attributes(graph, degrees, "degree")
    return clusters, degrees


def layout_nodes(
    nodes: list[str],
    edges: dict[tuple[str, str], dict[str, float | str]],
    width: int,
    height: int,
    seed: int,
) -> dict[str, tuple[float, float]]:
    rng = random.Random(seed)
    margin = 120
    positions = {
        node: (
            rng.uniform(margin, width - margin),
            rng.uniform(margin, height - margin),
        )
        for node in nodes
    }
    if len(nodes) <= 1:
        return positions

    area = (width - 2 * margin) * (height - 2 * margin)
    k = math.sqrt(area / len(nodes))
    temperature = min(width, height) / 8
    adjacency = set(edges.keys())

    for _ in range(250):
        disp = {node: [0.0, 0.0] for node in nodes}

        for i, v in enumerate(nodes):
            vx, vy = positions[v]
            for u in nodes[i + 1 :]:
                ux, uy = positions[u]
                dx = vx - ux
                dy = vy - uy
                dist = math.hypot(dx, dy) or 0.01
                force = (k * k) / dist
                fx = dx / dist * force
                fy = dy / dist * force
                disp[v][0] += fx
                disp[v][1] += fy
                disp[u][0] -= fx
                disp[u][1] -= fy

        for source, target in adjacency:
            sx, sy = positions[source]
            tx, ty = positions[target]
            dx = sx - tx
            dy = sy - ty
            dist = math.hypot(dx, dy) or 0.01
            force = (dist * dist) / k
            fx = dx / dist * force
            fy = dy / dist * force
            disp[source][0] -= fx
            disp[source][1] -= fy
            disp[target][0] += fx
            disp[target][1] += fy

        for node in nodes:
            dx, dy = disp[node]
            dist = math.hypot(dx, dy) or 0.01
            x, y = positions[node]
            x += dx / dist * min(dist, temperature)
            y += dy / dist * min(dist, temperature)
            positions[node] = (
                min(width - margin, max(margin, x)),
                min(height - margin, max(margin, y)),
            )
        temperature *= 0.97

    return scale_positions(positions, width, height, margin)


def scale_positions(
    positions: dict[str, tuple[float, float]], width: int, height: int, margin: int
) -> dict[str, tuple[float, float]]:
    if not positions:
        return positions
    xs = [value[0] for value in positions.values()]
    ys = [value[1] for value in positions.values()]
    min_x, max_x = min(xs), max(xs)
    min_y, max_y = min(ys), max(ys)
    span_x = max(max_x - min_x, 1.0)
    span_y = max(max_y - min_y, 1.0)
    scaled = {}
    for node, (x, y) in positions.items():
        scaled[node] = (
            margin + (x - min_x) / span_x * (width - 2 * margin),
            margin + (y - min_y) / span_y * (height - 2 * margin),
        )
    return scaled


def set_pixel(image: bytearray, width: int, height: int, x: int, y: int, color: tuple[int, int, int]) -> None:
    if 0 <= x < width and 0 <= y < height:
        offset = (y * width + x) * 3
        image[offset : offset + 3] = bytes(color)


def draw_line(
    image: bytearray,
    width: int,
    height: int,
    x0: float,
    y0: float,
    x1: float,
    y1: float,
    color: tuple[int, int, int],
) -> None:
    steps = max(1, int(math.hypot(x1 - x0, y1 - y0)))
    for step in range(steps + 1):
        t = step / steps
        x = int(round(x0 + (x1 - x0) * t))
        y = int(round(y0 + (y1 - y0) * t))
        set_pixel(image, width, height, x, y, color)


def draw_circle(
    image: bytearray,
    width: int,
    height: int,
    cx: float,
    cy: float,
    radius: int,
    color: tuple[int, int, int],
) -> None:
    x_center = int(round(cx))
    y_center = int(round(cy))
    r2 = radius * radius
    for y in range(y_center - radius, y_center + radius + 1):
        for x in range(x_center - radius, x_center + radius + 1):
            if (x - x_center) ** 2 + (y - y_center) ** 2 <= r2:
                set_pixel(image, width, height, x, y, color)


def write_png(path: Path, width: int, height: int, image: bytearray) -> None:
    raw = bytearray()
    for y in range(height):
        raw.append(0)
        start = y * width * 3
        raw.extend(image[start : start + width * 3])

    def chunk(name: bytes, data: bytes) -> bytes:
        payload = name + data
        return struct.pack(">I", len(data)) + payload + struct.pack(">I", zlib.crc32(payload) & 0xFFFFFFFF)

    png = b"\x89PNG\r\n\x1a\n"
    png += chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0))
    png += chunk(b"IDAT", zlib.compress(bytes(raw), level=9))
    png += chunk(b"IEND", b"")
    path.write_bytes(png)


def render_png(
    path: Path,
    nodes: list[str],
    edges: dict[tuple[str, str], dict[str, float | str]],
    positions: dict[str, tuple[float, float]],
    clusters: dict[str, int],
    degrees: dict[str, int],
    width: int,
    height: int,
) -> None:
    image = bytearray([255, 255, 255] * width * height)
    for source, target in edges:
        sx, sy = positions[source]
        tx, ty = positions[target]
        draw_line(image, width, height, sx, sy, tx, ty, (215, 219, 224))

    max_degree = max(degrees.values()) if degrees else 1
    for node in nodes:
        color = PALETTE[(clusters[node] - 1) % len(PALETTE)]
        radius = 8 + int(12 * math.sqrt(degrees[node] / max_degree)) if max_degree else 8
        x, y = positions[node]
        draw_circle(image, width, height, x, y, radius + 2, (255, 255, 255))
        draw_circle(image, width, height, x, y, radius, color)
    write_png(path, width, height, image)


def write_tables(
    outdir: Path,
    nodes: list[str],
    metadata: dict[str, dict[str, str]],
    edges: dict[tuple[str, str], dict[str, float | str]],
    clusters: dict[str, int],
    degrees: dict[str, int],
) -> None:
    node_path = outdir / "network_nodes.tsv"
    edge_path = outdir / "network_edges.tsv"
    with node_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "accession",
            "cluster",
            "degree",
            "protein_name",
            "organism",
            "length",
            "reviewed",
            "ec_number",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for node in nodes:
            row = metadata[node]
            writer.writerow(
                {
                    "accession": node,
                    "cluster": clusters[node],
                    "degree": degrees[node],
                    "protein_name": row.get("Protein names", ""),
                    "organism": row.get("Organism", ""),
                    "length": row.get("Length", ""),
                    "reviewed": row.get("Reviewed", ""),
                    "ec_number": row.get("EC number", ""),
                }
            )
    with edge_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = ["source", "target", "pident", "coverage", "evalue", "score", "bitscore"]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for edge in edges.values():
            writer.writerow(edge)


def write_graphml(
    path: Path,
    nodes: list[str],
    edges: dict[tuple[str, str], dict[str, float | str]],
    clusters: dict[str, int],
    degrees: dict[str, int],
) -> None:
    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">',
        '<key id="cluster" for="node" attr.name="cluster" attr.type="int"/>',
        '<key id="degree" for="node" attr.name="degree" attr.type="int"/>',
        '<key id="pident" for="edge" attr.name="pident" attr.type="double"/>',
        '<key id="coverage" for="edge" attr.name="coverage" attr.type="double"/>',
        '<key id="evalue" for="edge" attr.name="evalue" attr.type="double"/>',
        '<key id="bitscore" for="edge" attr.name="bitscore" attr.type="double"/>',
        '<graph id="SSN" edgedefault="undirected">',
    ]
    for node in nodes:
        lines.extend(
            [
                f'<node id="{node}">',
                f'<data key="cluster">{clusters[node]}</data>',
                f'<data key="degree">{degrees[node]}</data>',
                "</node>",
            ]
        )
    for index, edge in enumerate(edges.values(), start=1):
        lines.extend(
            [
                f'<edge id="e{index}" source="{edge["source"]}" target="{edge["target"]}">',
                f'<data key="pident">{edge["pident"]}</data>',
                f'<data key="coverage">{edge["coverage"]}</data>',
                f'<data key="evalue">{edge["evalue"]}</data>',
                f'<data key="bitscore">{edge["bitscore"]}</data>',
                "</edge>",
            ]
        )
    lines.extend(["</graph>", "</graphml>", ""])
    path.write_text("\n".join(lines), encoding="utf-8")


def dot_quote(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def html_escape(value: str) -> str:
    return (
        value.replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
        .replace('"', "&quot;")
    )


def cluster_sizes(clusters: dict[str, int]) -> dict[int, int]:
    sizes: dict[int, int] = defaultdict(int)
    for cluster in clusters.values():
        sizes[cluster] += 1
    return dict(sizes)


def cluster_color(cluster: int, featured_clusters: set[int]) -> str:
    if cluster not in featured_clusters:
        return "#b8bec6"
    return HEX_PALETTE[(cluster - 1) % len(HEX_PALETTE)]


def legend_html(
    sizes: dict[int, int],
    featured_clusters: list[int],
    min_identity: float,
    min_coverage: float,
    max_evalue: float,
) -> str:
    rows = [
        "<TR><TD ALIGN=\"LEFT\" COLSPAN=\"2\"><B>SSN legend</B></TD></TR>",
        f"<TR><TD ALIGN=\"LEFT\" COLSPAN=\"2\">Identity &gt;= {min_identity:g}%</TD></TR>",
        f"<TR><TD ALIGN=\"LEFT\" COLSPAN=\"2\">Coverage &gt;= {min_coverage:g}</TD></TR>",
        f"<TR><TD ALIGN=\"LEFT\" COLSPAN=\"2\">E-value &lt;= {max_evalue:g}</TD></TR>",
        "<TR><TD ALIGN=\"LEFT\" COLSPAN=\"2\">Node label: UniProt accession</TD></TR>",
        "<TR><TD ALIGN=\"LEFT\" COLSPAN=\"2\">Node size: degree</TD></TR>",
        "<TR><TD ALIGN=\"LEFT\" COLSPAN=\"2\">Node color: component</TD></TR>",
        "<TR><TD COLSPAN=\"2\"></TD></TR>",
    ]
    for cluster in featured_clusters:
        color = cluster_color(cluster, set(featured_clusters))
        rows.append(
            f"<TR><TD WIDTH=\"16\" BGCOLOR=\"{color}\"></TD><TD ALIGN=\"LEFT\">Component {cluster} (n={sizes[cluster]})</TD></TR>"
        )
    other_count = sum(size for cluster, size in sizes.items() if cluster not in featured_clusters)
    if other_count:
        rows.append(
            f"<TR><TD WIDTH=\"16\" BGCOLOR=\"#b8bec6\"></TD><TD ALIGN=\"LEFT\">Other components (n={other_count})</TD></TR>"
        )
    return "<<TABLE BORDER=\"0\" CELLBORDER=\"0\" CELLSPACING=\"2\">" + "".join(rows) + "</TABLE>>"


def write_dot(
    path: Path,
    nodes: list[str],
    edges: dict[tuple[str, str], dict[str, float | str]],
    clusters: dict[str, int],
    degrees: dict[str, int],
    min_identity: float,
    min_coverage: float,
    max_evalue: float,
    legend_components: int,
    label_font_size: int,
) -> None:
    max_degree = max(degrees.values()) if degrees else 1
    sizes = cluster_sizes(clusters)
    featured_clusters = sorted(sizes, key=lambda cluster: (-sizes[cluster], cluster))[:legend_components]
    featured_set = set(featured_clusters)
    lines = [
        "graph SSN {",
        "  graph [layout=sfdp, overlap=prism, pack=true, packmode=\"graph\", sep=\"+35\", K=1.4, repulsiveforce=2.2, splines=true, outputorder=edgesfirst, bgcolor=white, labelloc=t, labeljust=l, pad=0.4];",
        "  node [fontname=Helvetica];",
        "  edge [color=\"#aeb7c1\", penwidth=1.4];",
        f"  graph_legend [shape=plain, label={legend_html(sizes, featured_clusters, min_identity, min_coverage, max_evalue)}];",
    ]
    for node in nodes:
        cluster = clusters[node]
        degree = degrees[node]
        fill = cluster_color(cluster, featured_set)
        size = 0.32 + 0.45 * math.sqrt(degree / max_degree) if max_degree else 0.32
        lines.append(
            "  "
            + dot_quote(node)
            + f" [shape=circle, label=\"\", xlabel={dot_quote(node)}, fontsize={label_font_size}, style=filled, penwidth=1.1, color=white, fillcolor=\"{fill}\", width={size:.3f}, height={size:.3f}, tooltip={dot_quote(node)}];"
        )
    for edge in edges.values():
        lines.append(
            "  "
            + dot_quote(str(edge["source"]))
            + " -- "
            + dot_quote(str(edge["target"]))
            + f" [weight={float(edge['bitscore']):.3f}, tooltip={dot_quote(str(edge['pident']))}];"
        )
    lines.append("}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def render_with_graphviz(dot_path: Path, svg_path: Path) -> str | None:
    engine = shutil.which("sfdp") or shutil.which("dot")
    if not engine:
        return None
    subprocess.run(
        [engine, "-Tsvg", str(dot_path), "-o", str(svg_path)],
        check=True,
    )
    return Path(engine).name


def graphviz_positions(dot_path: Path, nodes: set[str]) -> dict[str, tuple[float, float]]:
    engine = shutil.which("sfdp") or shutil.which("dot")
    if not engine:
        return {}
    result = subprocess.run(
        [engine, "-Tplain", str(dot_path)],
        check=True,
        capture_output=True,
        text=True,
    )
    positions: dict[str, tuple[float, float]] = {}
    for line in result.stdout.splitlines():
        fields = line.split()
        if len(fields) >= 4 and fields[0] == "node" and fields[1] in nodes:
            positions[fields[1]] = (float(fields[2]) * 72.0, -float(fields[3]) * 72.0)
    return positions


def add_att(parent: ET.Element, name: str, value: object, value_type: str = "string") -> None:
    ET.SubElement(
        parent,
        "att",
        {
            "name": name,
            "type": value_type,
            "value": "" if value is None else str(value),
        },
    )


def infer_xgmml_type(value: object) -> str:
    if isinstance(value, bool):
        return "boolean"
    if isinstance(value, int):
        return "integer"
    if isinstance(value, float):
        return "real"
    return "string"


def write_xgmml(
    path: Path,
    nodes: list[str],
    metadata: dict[str, dict[str, str]],
    edges: dict[tuple[str, str], dict[str, float | str]],
    clusters: dict[str, int],
    degrees: dict[str, int],
    positions: dict[str, tuple[float, float]],
    legend_components: int,
) -> None:
    ET.register_namespace("", "http://www.cs.rpi.edu/XGMML")
    ET.register_namespace("cy", "http://www.cytoscape.org")
    graph = ET.Element(
        "graph",
        {
            "label": "Sequence Similarity Network",
            "directed": "0",
            "xmlns": "http://www.cs.rpi.edu/XGMML",
            "xmlns:cy": "http://www.cytoscape.org",
        },
    )
    add_att(graph, "description", "Protein sequence similarity network generated from all-vs-all BLASTP.")
    sizes = cluster_sizes(clusters)
    featured_clusters = sorted(sizes, key=lambda cluster: (-sizes[cluster], cluster))[:legend_components]
    featured_set = set(featured_clusters)
    max_degree = max(degrees.values()) if degrees else 1

    for node in nodes:
        row = metadata[node]
        cluster = clusters[node]
        degree = degrees[node]
        color = cluster_color(cluster, featured_set)
        size = 28.0 + 32.0 * math.sqrt(degree / max_degree) if max_degree else 28.0
        x, y = positions.get(node, (0.0, 0.0))
        node_el = ET.SubElement(graph, "node", {"id": node, "label": node})
        node_attrs = {
            "accession": node,
            "protein_name": row.get("Protein names", ""),
            "organism": row.get("Organism", ""),
            "length": int(row["Length"]) if row.get("Length", "").isdigit() else row.get("Length", ""),
            "reviewed": row.get("Reviewed", ""),
            "ec_number": row.get("EC number", ""),
            "cluster": cluster,
            "degree": degree,
            "component_size": sizes[cluster],
        }
        for name, value in node_attrs.items():
            add_att(node_el, name, value, infer_xgmml_type(value))
        ET.SubElement(
            node_el,
            "graphics",
            {
                "type": "ELLIPSE",
                "x": f"{x:.3f}",
                "y": f"{y:.3f}",
                "w": f"{size:.3f}",
                "h": f"{size:.3f}",
                "fill": color,
                "outline": "#ffffff",
                "width": "1.0",
            },
        )

    for index, edge in enumerate(edges.values(), start=1):
        source = str(edge["source"])
        target = str(edge["target"])
        edge_el = ET.SubElement(
            graph,
            "edge",
            {
                "id": f"e{index}",
                "label": f"{source} -- {target}",
                "source": source,
                "target": target,
            },
        )
        for name in ["pident", "coverage", "evalue", "score", "bitscore"]:
            add_att(edge_el, name, edge[name], "real")
        ET.SubElement(edge_el, "graphics", {"fill": "#aeb7c1", "width": "1.4"})

    tree = ET.ElementTree(graph)
    ET.indent(tree, space="  ")
    tree.write(path, encoding="utf-8", xml_declaration=True)


def main() -> int:
    args = parse_args()
    metadata = read_metadata(args.metadata)
    nodes = sorted(metadata)
    edges = read_edges(args.similarity, args.min_identity, args.min_coverage, args.max_evalue)
    graph = build_networkx_graph(nodes, metadata, edges)
    if graph is not None:
        clusters, degrees = networkx_components_and_degrees(graph)
        graph_logic = f"NetworkX {nx.__version__}"
    else:
        clusters = connected_components(nodes, edges)
        degrees = degree_by_node(nodes, edges)
        graph_logic = "built-in fallback"

    args.outdir.mkdir(parents=True, exist_ok=True)
    write_tables(args.outdir, nodes, metadata, edges, clusters, degrees)
    if graph is not None:
        nx.write_graphml(graph, args.outdir / "network.graphml")
    else:
        write_graphml(args.outdir / "network.graphml", nodes, edges, clusters, degrees)
    dot_path = args.outdir / "network.dot"
    png_path = args.outdir / "network.png"
    svg_path = args.outdir / "network.svg"
    write_dot(
        dot_path,
        nodes,
        edges,
        clusters,
        degrees,
        args.min_identity,
        args.min_coverage,
        args.max_evalue,
        args.legend_components,
        args.label_font_size,
    )
    graphviz_engine = render_with_graphviz(dot_path, svg_path)
    if graphviz_engine is None:
        positions = layout_nodes(nodes, edges, args.width, args.height, args.seed)
        render_png(png_path, nodes, edges, positions, clusters, degrees, args.width, args.height)
        rendering = "built-in fallback PNG renderer"
    else:
        positions = graphviz_positions(dot_path, set(nodes))
        rendering = f"Graphviz {graphviz_engine}"
    xgmml_path = args.outdir / "network.xgmml"
    write_xgmml(
        xgmml_path,
        nodes,
        metadata,
        edges,
        clusters,
        degrees,
        positions,
        args.legend_components,
    )

    component_sizes = cluster_sizes(clusters)
    summary = "\n".join(
        [
            "Sequence similarity network summary",
            "===================================",
            f"Similarity table: {args.similarity}",
            f"Metadata table: {args.metadata}",
            f"Minimum identity: {args.min_identity}",
            f"Minimum coverage: {args.min_coverage}",
            f"Maximum e-value: {args.max_evalue}",
            f"Nodes: {len(nodes)}",
            f"Edges: {len(edges)}",
            f"Connected components: {len(component_sizes)}",
            f"Graph logic: {graph_logic}",
            f"Rendering: {rendering}",
            "Component sizes: "
            + ", ".join(f"{cluster}:{size}" for cluster, size in sorted(component_sizes.items())),
            f"SVG figure: {svg_path if graphviz_engine is not None else 'not generated; Graphviz unavailable'}",
            f"PNG fallback figure: {png_path if graphviz_engine is None else 'not generated; Graphviz SVG render available'}",
            f"GraphML: {args.outdir / 'network.graphml'}",
            f"DOT: {dot_path}",
            f"XGMML: {xgmml_path}",
            "",
        ]
    )
    (args.outdir / "network_summary.txt").write_text(summary, encoding="utf-8")
    print(summary)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)
