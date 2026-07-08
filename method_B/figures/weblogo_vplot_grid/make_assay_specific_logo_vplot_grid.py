#!/usr/bin/env python3
"""Assemble assay-specific TF weblogos above matched V-plot panels."""
from __future__ import annotations

import csv
from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent
SCATTER_ROOT = ROOT / "v_apex_scores" / "methodB_scatter_best_assay_specific"
PANEL_ROOT = SCATTER_ROOT / "curated_vplot_panels"
OUTDIR = PANEL_ROOT / "weblogo_top_vplot_bottom"

LOGO_DIRS = [
    ROOT / "weblogo_out" / "logo_scatter_present_screen" / "logo",
    ROOT / "weblogo_out" / "logo_diverse_screen" / "logo",
    ROOT / "weblogo_out" / "logo_verified_strand_aligned" / "logo",
    ROOT / "weblogo_out" / "logo_verified_NFE2_aligned" / "logo",
    ROOT / "weblogo_out" / "logo_curated_AT_balanced" / "logo",
    ROOT / "weblogo_out" / "logo_from_motif_coords" / "logo",
    ROOT / "weblogo_out" / "logo" / "",
]

V_PANELS = [
    PANEL_ROOT / "diverse_selected_panels",
    PANEL_ROOT / "from_server_channel_panels",
    PANEL_ROOT / "replacement_candidate_panels",
]

ASSAYS = {
    "loMNase": ["TFE3", "NFIC"],
    "DNase": ["REST", "NFE2"],
    "ATAC": ["CEBPG", "CREB1"],
}

AT_GC = {
    "NFIC": ("0.4091", "0.5909"),
    "TFE3": ("0.4825", "0.5175"),
    "MAFG": ("0.4958", "0.5042"),
    "NFE2": ("0.5307", "0.4693"),
    "CEBPG": ("0.6328", "0.3672"),
    "MLX": ("0.4129", "0.5871"),
    "CREB1": ("0.6092", "0.3908"),
    "REST": ("0.3828", "0.6172"),
    "JUND": ("0.5615", "0.4385"),
    "MECOM": ("0.7433", "0.2567"),
    "TBP": ("0.7362", "0.2638"),
    "RFX1": ("0.4737", "0.5263"),
    "ZNF324": ("0.4854", "0.5146"),
}

CANVAS_BG = "white"
COL_W = 760
COL_GAP = 90
LEFT_MARGIN = 155
RIGHT_MARGIN = 55
TOP_MARGIN = 42
BOTTOM_MARGIN = 130
ASSAY_TITLE_H = 56
TF_TITLE_H = 44
LOGO_H = 185
LOGO_VPLOT_GAP = 18
VPLOT_W = 650
VPLOT_BOX_H = 545
PAIR_GAP = 86

LOGO_TITLE_CROP_TOP = 56
VPLOT_TITLE_CROP_FALLBACK = 78


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    candidates = [
        Path("C:/Windows/Fonts/arialbd.ttf") if bold else Path("C:/Windows/Fonts/arial.ttf"),
        Path("C:/Windows/Fonts/calibrib.ttf") if bold else Path("C:/Windows/Fonts/calibri.ttf"),
        Path("C:/Windows/Fonts/arial.ttf"),
    ]
    for path in candidates:
        if path.exists():
            return ImageFont.truetype(str(path), size)
    return ImageFont.load_default()


def resize_by_height(img: Image.Image, height: int) -> Image.Image:
    width = round(img.width * height / img.height)
    return img.resize((width, height), Image.Resampling.LANCZOS)


def resize_by_width(img: Image.Image, width: int) -> Image.Image:
    height = round(img.height * width / img.width)
    return img.resize((width, height), Image.Resampling.LANCZOS)


def find_logo(tf: str) -> Path:
    for directory in LOGO_DIRS:
        path = directory / f"{tf}.png"
        if path.exists():
            return path
    raise FileNotFoundError(f"Missing logo for {tf}")


def find_vplot(assay: str, tf: str) -> Path:
    for root in V_PANELS:
        path = root / assay / f"{tf}.png"
        if path.exists():
            return path
    raise FileNotFoundError(f"Missing V-plot for {assay} {tf}")


def prepare_logo(tf: str) -> Image.Image:
    img = Image.open(find_logo(tf)).convert("RGB")
    img = img.crop((0, LOGO_TITLE_CROP_TOP, img.width, img.height))
    return resize_by_height(img, LOGO_H)


def prepare_vplot(assay: str, tf: str) -> Image.Image:
    img = Image.open(find_vplot(assay, tf)).convert("RGB")
    img = img.crop((0, detect_vplot_crop_top(img), img.width, img.height))
    return resize_by_width(img, VPLOT_W)


def detect_vplot_crop_top(img: Image.Image) -> int:
    """Crop just above the plot frame so baked-in panel titles disappear."""
    pix = img.load()
    threshold = int(img.width * 0.35)
    for y in range(min(180, img.height)):
        dark = 0
        for x in range(img.width):
            r, g, b = pix[x, y]
            if r < 80 and g < 80 and b < 80:
                dark += 1
        if dark > threshold:
            return max(0, y - 24)
    return VPLOT_TITLE_CROP_FALLBACK


def draw_centered_text(
    draw: ImageDraw.ImageDraw,
    xy: tuple[int, int, int, int],
    text: str,
    fill: str,
    font_obj: ImageFont.ImageFont,
) -> None:
    x0, y0, x1, y1 = xy
    bb = draw.textbbox((0, 0), text, font=font_obj)
    w = bb[2] - bb[0]
    h = bb[3] - bb[1]
    draw.text((x0 + (x1 - x0 - w) / 2, y0 + (y1 - y0 - h) / 2 - bb[1]), text, fill=fill, font=font_obj)


def paste_centered(canvas: Image.Image, img: Image.Image, box: tuple[int, int, int, int]) -> None:
    x0, y0, x1, y1 = box
    x = round(x0 + (x1 - x0 - img.width) / 2)
    y = round(y0 + (y1 - y0 - img.height) / 2)
    canvas.paste(img, (x, y))


def main() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)

    assay_font = font(42, bold=True)
    tf_font = font(32)
    label_font = font(34)

    pair_h = TF_TITLE_H + LOGO_H + LOGO_VPLOT_GAP + VPLOT_BOX_H
    width = LEFT_MARGIN + len(ASSAYS) * COL_W + (len(ASSAYS) - 1) * COL_GAP + RIGHT_MARGIN
    height = TOP_MARGIN + ASSAY_TITLE_H + 2 * pair_h + PAIR_GAP + BOTTOM_MARGIN

    canvas = Image.new("RGB", (width, height), CANVAS_BG)
    draw = ImageDraw.Draw(canvas)

    manifest_rows = []
    for col, (assay, tfs) in enumerate(ASSAYS.items()):
        col_x = LEFT_MARGIN + col * (COL_W + COL_GAP)
        draw_centered_text(
            draw,
            (col_x, TOP_MARGIN - 4, col_x + COL_W, TOP_MARGIN + ASSAY_TITLE_H),
            assay,
            "black",
            assay_font,
        )
        y = TOP_MARGIN + ASSAY_TITLE_H
        for row, tf in enumerate(tfs):
            pair_y = y + row * (pair_h + PAIR_GAP)
            draw_centered_text(
                draw,
                (col_x, pair_y, col_x + COL_W, pair_y + TF_TITLE_H),
                tf,
                "black",
                tf_font,
            )

            logo = prepare_logo(tf)
            logo_box = (
                col_x,
                pair_y + TF_TITLE_H,
                col_x + COL_W,
                pair_y + TF_TITLE_H + LOGO_H,
            )
            paste_centered(canvas, logo, logo_box)

            vplot = prepare_vplot(assay, tf)
            vplot_box = (
                col_x,
                pair_y + TF_TITLE_H + LOGO_H + LOGO_VPLOT_GAP,
                col_x + COL_W,
                pair_y + TF_TITLE_H + LOGO_H + LOGO_VPLOT_GAP + VPLOT_BOX_H,
            )
            paste_centered(canvas, vplot, vplot_box)

            at, gc = AT_GC.get(tf, ("", ""))
            manifest_rows.append(
                {
                    "assay": assay,
                    "TF": tf,
                    "logo_png": str(find_logo(tf)),
                    "vplot_png": str(find_vplot(assay, tf)),
                    "AT": at,
                    "GC": gc,
                }
            )

    x_label = "Distance from motif (bp)"
    y_label = "Fragment length (bp)"
    draw_centered_text(draw, (0, height - 82, width, height - 28), x_label, "black", label_font)

    tmp = Image.new("RGBA", (height, 80), (255, 255, 255, 0))
    tmp_draw = ImageDraw.Draw(tmp)
    bb = tmp_draw.textbbox((0, 0), y_label, font=label_font)
    tmp_draw.text(((height - (bb[2] - bb[0])) / 2, 8), y_label, fill="black", font=label_font)
    rotated = tmp.rotate(90, expand=True)
    canvas.paste(rotated.convert("RGB"), (24, (height - rotated.height) // 2), rotated)

    out_png = OUTDIR / "assay_specific_weblogo_above_vplot.png"
    out_pdf = OUTDIR / "assay_specific_weblogo_above_vplot.pdf"
    canvas.save(out_png, optimize=True)
    canvas.save(out_pdf, "PDF", resolution=150)

    with (OUTDIR / "assay_specific_weblogo_above_vplot_manifest.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["assay", "TF", "logo_png", "vplot_png", "AT", "GC"])
        writer.writeheader()
        writer.writerows(manifest_rows)

    print(f"wrote {out_png}")
    print(f"wrote {out_pdf}")


if __name__ == "__main__":
    main()
