import re
import argparse
from io import BytesIO
from pathlib import Path

from pypdf import PdfReader, PdfWriter
from pypdf.generic import DecodedStreamObject, NameObject
from reportlab.pdfgen import canvas
from reportlab.lib.colors import HexColor, white


PROJECT_ROOT = Path(__file__).resolve().parents[4]
DEFAULT_INPUT_PDF = Path("/Users/liuz/Downloads/95000bc2-c6e7-4767-af7d-37c76e7a2175.pdf")
DEFAULT_OUTPUT_PDF = (
    PROJECT_ROOT
    / "final"
    / "figure_7"
    / "mp7_myeloid_spatial_recolored"
    / "outputs"
    / "figure7_spatial_mp7_myeloid_recolored.pdf"
)

TARGET_COLORS = {
    "MP_7": "#C52B2F",
    "Mono_CDC27": "#0072B2",
    "Mono_SLC2A3": "#009E73",
    "Neutrophil": "#E69F00",
    "other cell types": "#BDBDBD",
}

# Source RGB triples from the uploaded Matplotlib categorical PDF.
# All non-highlight categories are intentionally collapsed to a shared gray.
SOURCE_TO_TARGET = {
    "1 1 0": TARGET_COLORS["other cell types"],  # A2ML1+ epi
    "0.1098039216 0.9019607843 1": TARGET_COLORS["other cell types"],  # B/Plasma
    "1 0.2039215686 1": TARGET_COLORS["other cell types"],  # CD4_CXCL13
    "1 0.2901960784 0.2745098039": TARGET_COLORS["other cell types"],  # CD4_Treg_CCR8
    "0 0.537254902 0.2549019608": TARGET_COLORS["other cell types"],  # CD4_Treg_FOXP3
    "0 0.4352941176 0.6509803922": TARGET_COLORS["other cell types"],  # CD8_Teff
    "0.6392156863 0 0.3490196078": TARGET_COLORS["other cell types"],  # CD8_Tex_PDCD1
    "1 0.8588235294 0.8980392157": TARGET_COLORS["other cell types"],  # CD8_prolif
    "0.4784313725 0.2862745098 0": TARGET_COLORS["other cell types"],  # MP_1
    "0 0 0.6509803922": TARGET_COLORS["other cell types"],  # MP_4
    "0.3882352941 1 0.6745098039": TARGET_COLORS["other cell types"],  # MP_5
    "0.7176470588 0.5921568627 0.3843137255": TARGET_COLORS["other cell types"],  # MP_6
    "0 0.3019607843 0.262745098": TARGET_COLORS["MP_7"],
    "0.5607843137 0.6901960784 1": TARGET_COLORS["other cell types"],  # Macro_CXCL5
    "0.6 0.4901960784 0.5294117647": TARGET_COLORS["other cell types"],  # Macro_CXCL9
    "0.3529411765 0 0.0274509804": TARGET_COLORS["Mono_CDC27"],
    "0.5019607843 0.5882352941 0.5764705882": TARGET_COLORS["Mono_SLC2A3"],
    "0.4156862745 0.2274509804 0.2980392157": TARGET_COLORS["Neutrophil"],
    "0.1058823529 0.2666666667 0": TARGET_COLORS["other cell types"],  # T/NK
    "0.3098039216 0.7764705882 0.0039215686": TARGET_COLORS["other cell types"],  # endo_PLVAP
    "0.231372549 0.3647058824 1": TARGET_COLORS["other cell types"],  # iCAF_CXCL12
    "0.2901960784 0.231372549 0.3254901961": TARGET_COLORS["other cell types"],  # iCAF_CXCL5
    "1 0.1843137255 0.5019607843": TARGET_COLORS["other cell types"],  # iCAF_CXCL6
    "0.3803921569 0.3803921569 0.3529411765": TARGET_COLORS["other cell types"],  # iCAF_TNC
    "0.7294117647 0.0352941176 0": TARGET_COLORS["other cell types"],  # myCAF_COL10A1
    "0.4196078431 0.4745098039 0": TARGET_COLORS["other cell types"],  # myCAF_DUX4
    "0 0.7607843137 0.6274509804": TARGET_COLORS["other cell types"],  # myCAF_MMP11
    "1 0.6666666667 0.5725490196": TARGET_COLORS["other cell types"],  # myCAF_TNFRSF21
    "1 0.5647058824 0.7882352941": TARGET_COLORS["other cell types"],  # vCAF
    "0.8666666667 0.937254902 1": TARGET_COLORS["other cell types"],
    "0.7254901961 0.0117647059 0.6666666667": TARGET_COLORS["other cell types"],
    "0.8196078431 0.3803921569 0": TARGET_COLORS["other cell types"],
    # Source colors from a second Matplotlib categorical palette used by related PDFs.
    "0.0078431373 0.2470588235 0.6470588235": TARGET_COLORS["other cell types"],  # A2ML1+ epi
    "0.4901960784 0.5294117647 0.7254901961": TARGET_COLORS["other cell types"],  # B/Plasma
    "0.7450980392 0.7568627451 0.831372549": TARGET_COLORS["other cell types"],  # CD4_CXCL13
    "0.8392156863 0.737254902 0.7529411765": TARGET_COLORS["other cell types"],  # CD4_Treg_CCR8
    "0.7333333333 0.4666666667 0.5176470588": TARGET_COLORS["other cell types"],  # CD8_Teff
    "0.5568627451 0.0235294118 0.231372549": TARGET_COLORS["other cell types"],  # CD8_prolif
    "0.2901960784 0.4352941176 0.8901960784": TARGET_COLORS["other cell types"],  # MP_1
    "0.5215686275 0.5843137255 0.8823529412": TARGET_COLORS["other cell types"],  # MP_2
    "0.7098039216 0.7333333333 0.8901960784": TARGET_COLORS["other cell types"],  # MP_3
    "0.9019607843 0.6862745098 0.7254901961": TARGET_COLORS["other cell types"],  # MP_4
    "0.8784313725 0.4823529412 0.568627451": TARGET_COLORS["other cell types"],  # MP_5
    "0.8274509804 0.2470588235 0.4156862745": TARGET_COLORS["other cell types"],  # MP_6
    "0.0666666667 0.7764705882 0.2196078431": TARGET_COLORS["MP_7"],
    "0.5529411765 0.8352941176 0.5764705882": TARGET_COLORS["other cell types"],  # Macro_CXCL5
    "0.7764705882 0.8705882353 0.7803921569": TARGET_COLORS["other cell types"],  # Macro_CXCL9
    "0.9176470588 0.8274509804 0.7764705882": TARGET_COLORS["Mono_CDC27"],
    "0.9411764706 0.7254901961 0.5529411765": TARGET_COLORS["Mono_SLC2A3"],
    "0.937254902 0.5921568627 0.031372549": TARGET_COLORS["Neutrophil"],
    "0.0588235294 0.8117647059 0.7529411765": TARGET_COLORS["other cell types"],  # T/NK
    "0.6117647059 0.8705882353 0.8392156863": TARGET_COLORS["other cell types"],  # cDC3_LAMP3
    "0.8352941176 0.9176470588 0.9058823529": TARGET_COLORS["other cell types"],  # iCAF_CXCL5
}


def hex_to_pdf_rgb(hex_color):
    hex_color = hex_color.lstrip("#")
    values = [int(hex_color[i : i + 2], 16) / 255 for i in (0, 2, 4)]
    return " ".join(f"{v:.10f}".rstrip("0").rstrip(".") for v in values)


def replace_page_content(page):
    contents = page.get_contents()
    if contents is None:
        return {}

    data = contents.get_data().decode("latin1")
    counts = {}
    for source, target_hex in SOURCE_TO_TARGET.items():
        target = hex_to_pdf_rgb(target_hex)
        pattern = r"(?<![\d.])" + r"\s+".join(re.escape(part) for part in source.split()) + r"(?![\d.])"
        data, counts[source] = re.subn(pattern, target, data)

    new_stream = DecodedStreamObject()
    new_stream.set_data(data.encode("latin1"))
    page[NameObject("/Contents")] = new_stream
    return counts


def make_legend_overlay(width, height, cover_x):
    packet = BytesIO()
    c = canvas.Canvas(packet, pagesize=(width, height))

    # Cover the original long legend, leaving the plotting panel untouched.
    c.setFillColor(white)
    c.setStrokeColor(white)
    c.rect(cover_x, 42, width - cover_x, height - 72, fill=1, stroke=0)
    c.setStrokeColorRGB(0, 0, 0)
    c.setLineWidth(0.5)
    c.line(cover_x, 42, cover_x, height - 38)

    legend_items = [
        ("MP_7", TARGET_COLORS["MP_7"]),
        ("Mono_CDC27", TARGET_COLORS["Mono_CDC27"]),
        ("Mono_SLC2A3", TARGET_COLORS["Mono_SLC2A3"]),
        ("Neutrophil", TARGET_COLORS["Neutrophil"]),
        ("other cell types", TARGET_COLORS["other cell types"]),
    ]

    x_dot = cover_x + 17
    x_text = cover_x + 30
    y = height - 98
    for label, color in legend_items:
        c.setFillColor(HexColor(color))
        c.circle(x_dot, y + 1.5, 3.5, fill=1, stroke=0)
        c.setFillColorRGB(0, 0, 0)
        c.setFont("Helvetica", 7.5)
        c.drawString(x_text, y - 1.5, label)
        y -= 14

    c.save()
    packet.seek(0)
    return PdfReader(packet).pages[0]


parser = argparse.ArgumentParser(description="Recolor MP7-myeloid spatial PDFs to the shared Figure 7 palette.")
parser.add_argument("--input-pdf", type=Path, default=DEFAULT_INPUT_PDF)
parser.add_argument("--output-pdf", type=Path, default=DEFAULT_OUTPUT_PDF)
parser.add_argument("--legend-cover-x", type=float, default=418)
args = parser.parse_args()

reader = PdfReader(args.input_pdf)
writer = PdfWriter()

all_counts = []
for page in reader.pages:
    counts = replace_page_content(page)
    all_counts.append(counts)
    width = float(page.mediabox.width)
    height = float(page.mediabox.height)
    page.merge_page(make_legend_overlay(width, height, args.legend_cover_x))
    writer.add_page(page)

args.output_pdf.parent.mkdir(parents=True, exist_ok=True)
with args.output_pdf.open("wb") as f:
    writer.write(f)

print(f"Saved: {args.output_pdf}")
print("Target colors:")
for label, color in TARGET_COLORS.items():
    print(f"  {label}: {color}")
print("Replacement counts:")
for i, counts in enumerate(all_counts, start=1):
    found = sum(counts.values())
    missed = [source for source, n in counts.items() if n == 0]
    print(f"  Page {i}: {found} replacements; {len(missed)} source colors missed")
    if missed:
        print("    missed:", missed)
