from __future__ import annotations

from pathlib import Path

from docx import Document
from docx.enum.section import WD_SECTION
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Pt, RGBColor


HERE = Path(__file__).resolve().parent
PROJECT = HERE.parent
EQUATION_DIR = HERE / "_forward_model_equations"
OUTPUT = HERE / "L波段前向亮温计算模型公式与文献来源.docx"
# Searchable Unicode equations are used instead of raster formulas.  This keeps
# the document self-contained on machines where the optional MathJax package is
# not installed, while Word still renders the equations with Cambria Math.
FORMULAS = {
    "model_input": "T_B = T_B(h_i, h_s, T_sfc, ice type, T_water, S_water)",
    "conductivity": "k_ice = 2.034 + 0.13 S_ice/(T_ice − 273.15)    ;    k_snow = 0.31 W·m⁻¹·K⁻¹",
    "heat_profile": "R_th = h_s/k_s + h_i/k_i    ;    q = (T_w − T_sfc)/R_th\nT(z) = T₀ + qz/k",
    "fyi_salinity": "S_FYI = 6.08 exp(−5.81h_i) + 7.409 exp(−0.5228h_i)    ;    S_FYI ≥ 1.5 ppt",
    "myi_salinity": "z_k = (k − 1/2)/N    ;    S_i,k = (S_max/2){1 − cos[π z_k^(a/(z_k+b))]}",
    "myi_constants": "S_max = 3.2 ppt    ;    a = 0.407    ;    b = 0.573",
    "snow_eps": "ε′_s = 1 + 0.7ρ_s + 0.7ρ_s²\nε″_s = 1.59×10⁶(0.52ρ_s + 0.62ρ_s²)(f⁻¹ + 1.23×10⁻¹⁴√f) exp(0.036T_s)\nε_s = ε′_s + iε″_s",
    "ice_eps": "ε_i = a₁ + a₂V_b + i(a₃ + a₄V_b)",
    "brine": "ρ_i = 0.917 − 1.403×10⁻⁴T_i    ;    V_b = 1000ρ_iS_i/[F₁(T_i) − ρ_iS_iF₂(T_i)]",
    "water_eps": "ε_w = ε_∞ + [ε_s(T,S) − ε_∞]/[1 + iωτ(T,S)] − iσ(T,S)/(ωε₀)    ;    ε_∞ = 4.9",
    "fresnel": "q_j = √(ε_j − sin²θ₀)\nr_H = (q₁ − q₂)/(q₁ + q₂)    ;    r_V = (ε₂q₁ − ε₁q₂)/(ε₂q₁ + ε₁q₂)    ;    R = |r|²",
    "attenuation": "τ_j = exp{−2k₀|Im[√(ε_j − sin²θ₀)]|d_j}    ;    k₀ = 2πf/c₀",
    "recursion": "A_j = (1 − τ_j)(1 + R_j,b τ_j)\nU_j = T_jA_j + (1 − R_j,b)U_(j+1) τ_j\nE_j = A_j + (1 − R_j,b)E_(j+1) τ_j",
    "top_surface": "T_B,surf = (1 − R_top)U₁    ;    e_eff = (1 − R_top)E₁",
    "toa": "T_B,TOA = T_B,surf + (1 − e_eff)T_B,cosm + ΔT_B,atm + ΔT_B,offset",
    "angular_mean": "T̄_B = [1/(2N_θ)] Σ_(p∈{H,V}) Σ_(θ=0°…40°) T_B(p,θ)",
}


def set_cell_shading(cell, fill: str) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    shd = tc_pr.find(qn("w:shd"))
    if shd is None:
        shd = OxmlElement("w:shd")
        tc_pr.append(shd)
    shd.set(qn("w:fill"), fill)


def set_cell_margins(cell, top=100, start=110, bottom=100, end=110) -> None:
    tc = cell._tc
    tc_pr = tc.get_or_add_tcPr()
    tc_mar = tc_pr.first_child_found_in("w:tcMar")
    if tc_mar is None:
        tc_mar = OxmlElement("w:tcMar")
        tc_pr.append(tc_mar)
    for margin, value in (("top", top), ("start", start), ("bottom", bottom), ("end", end)):
        node = tc_mar.find(qn(f"w:{margin}"))
        if node is None:
            node = OxmlElement(f"w:{margin}")
            tc_mar.append(node)
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")


def set_table_borders(table, color="D9D9D9", size="6") -> None:
    tbl_pr = table._tbl.tblPr
    borders = tbl_pr.first_child_found_in("w:tblBorders")
    if borders is None:
        borders = OxmlElement("w:tblBorders")
        tbl_pr.append(borders)
    for edge in ("top", "left", "bottom", "right", "insideH", "insideV"):
        element = borders.find(qn(f"w:{edge}"))
        if element is None:
            element = OxmlElement(f"w:{edge}")
            borders.append(element)
        element.set(qn("w:val"), "single")
        element.set(qn("w:sz"), size)
        element.set(qn("w:space"), "0")
        element.set(qn("w:color"), color)


def set_repeat_table_header(row) -> None:
    tr_pr = row._tr.get_or_add_trPr()
    repeat = OxmlElement("w:tblHeader")
    repeat.set(qn("w:val"), "true")
    tr_pr.append(repeat)


def set_run_font(run, latin="Times New Roman", east_asia="宋体", size=10.5, bold=None) -> None:
    run.font.name = latin
    run._element.get_or_add_rPr().rFonts.set(qn("w:ascii"), latin)
    run._element.get_or_add_rPr().rFonts.set(qn("w:hAnsi"), latin)
    run._element.get_or_add_rPr().rFonts.set(qn("w:eastAsia"), east_asia)
    run.font.size = Pt(size)
    if bold is not None:
        run.bold = bold


def add_hyperlink(paragraph, text: str, url: str):
    part = paragraph.part
    rel_id = part.relate_to(url, "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink", is_external=True)
    hyperlink = OxmlElement("w:hyperlink")
    hyperlink.set(qn("r:id"), rel_id)
    run = OxmlElement("w:r")
    run_pr = OxmlElement("w:rPr")
    color = OxmlElement("w:color")
    color.set(qn("w:val"), "0563C1")
    underline = OxmlElement("w:u")
    underline.set(qn("w:val"), "single")
    rfonts = OxmlElement("w:rFonts")
    rfonts.set(qn("w:ascii"), "Times New Roman")
    rfonts.set(qn("w:hAnsi"), "Times New Roman")
    rfonts.set(qn("w:eastAsia"), "宋体")
    run_pr.extend([rfonts, color, underline])
    run.append(run_pr)
    text_node = OxmlElement("w:t")
    text_node.text = text
    run.append(text_node)
    hyperlink.append(run)
    paragraph._p.append(hyperlink)
    return hyperlink


def set_keep_with_next(paragraph, value=True) -> None:
    paragraph.paragraph_format.keep_with_next = value


def add_body(doc, text: str, bold_lead: str | None = None, after=5):
    p = doc.add_paragraph(style="Body Text")
    p.paragraph_format.space_after = Pt(after)
    p.paragraph_format.line_spacing = 1.25
    if bold_lead and text.startswith(bold_lead):
        lead = p.add_run(bold_lead)
        set_run_font(lead, bold=True)
        run = p.add_run(text[len(bold_lead):])
        set_run_font(run)
    else:
        run = p.add_run(text)
        set_run_font(run)
    return p


def add_source(doc, text: str):
    p = doc.add_paragraph()
    p.paragraph_format.left_indent = Inches(0.22)
    p.paragraph_format.space_after = Pt(7)
    p.paragraph_format.keep_together = True
    lead = p.add_run("来源定位  ")
    set_run_font(lead, size=9, bold=True)
    lead.font.color.rgb = RGBColor(64, 64, 64)
    run = p.add_run(text)
    set_run_font(run, size=9)
    run.font.color.rgb = RGBColor(64, 64, 64)
    return p


def add_formula(doc, caption: str, image_name: str, width: float, source: str):
    caption_p = doc.add_paragraph()
    caption_p.paragraph_format.space_before = Pt(3)
    caption_p.paragraph_format.space_after = Pt(2)
    set_keep_with_next(caption_p)
    run = caption_p.add_run(caption)
    set_run_font(run, size=9.5, bold=True)
    run.font.color.rgb = RGBColor(31, 78, 121)

    p = doc.add_paragraph()
    p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    p.paragraph_format.space_after = Pt(2)
    p.paragraph_format.keep_together = True
    p.paragraph_format.left_indent = Inches(0.12)
    p.paragraph_format.right_indent = Inches(0.12)
    p.paragraph_format.line_spacing = 1.12
    p_pr = p._p.get_or_add_pPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:fill"), "F7F9FB")
    p_pr.append(shd)
    run = p.add_run(FORMULAS[image_name])
    set_run_font(run, latin="Cambria Math", east_asia="Cambria Math", size=10.5)
    add_source(doc, source)


def add_heading(doc, text: str, level: int):
    p = doc.add_heading(text, level=level)
    p.paragraph_format.space_before = Pt(10 if level == 1 else 7)
    p.paragraph_format.space_after = Pt(5)
    set_keep_with_next(p)
    return p


def add_bullet(doc, text: str, level=0):
    style = "List Bullet" if level == 0 else "List Bullet 2"
    p = doc.add_paragraph(style=style)
    p.paragraph_format.space_after = Pt(3)
    p.paragraph_format.line_spacing = 1.15
    set_run_font(p.add_run(text))
    return p


def set_col_width(cell, width_inches: float) -> None:
    tc_pr = cell._tc.get_or_add_tcPr()
    tc_w = tc_pr.find(qn("w:tcW"))
    if tc_w is None:
        tc_w = OxmlElement("w:tcW")
        tc_pr.append(tc_w)
    tc_w.set(qn("w:w"), str(int(width_inches * 1440)))
    tc_w.set(qn("w:type"), "dxa")


def populate_table(table, headers, rows, widths, font_size=8.6):
    table.autofit = False
    set_table_borders(table)
    header = table.rows[0]
    set_repeat_table_header(header)
    for idx, text in enumerate(headers):
        cell = header.cells[idx]
        cell.text = text
        set_cell_shading(cell, "1F4E79")
        cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        set_cell_margins(cell)
        set_col_width(cell, widths[idx])
        for paragraph in cell.paragraphs:
            paragraph.alignment = WD_ALIGN_PARAGRAPH.CENTER
            paragraph.paragraph_format.space_after = Pt(0)
            for run in paragraph.runs:
                set_run_font(run, size=font_size, bold=True)
                run.font.color.rgb = RGBColor(255, 255, 255)
    for row_idx, values in enumerate(rows, start=1):
        row = table.rows[row_idx]
        for col_idx, value in enumerate(values):
            cell = row.cells[col_idx]
            cell.text = str(value)
            if row_idx % 2 == 0:
                set_cell_shading(cell, "F3F6F9")
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            set_cell_margins(cell)
            set_col_width(cell, widths[col_idx])
            for paragraph in cell.paragraphs:
                paragraph.paragraph_format.space_after = Pt(0)
                paragraph.paragraph_format.line_spacing = 1.05
                paragraph.alignment = WD_ALIGN_PARAGRAPH.LEFT if col_idx in (1, 2, 3) else WD_ALIGN_PARAGRAPH.CENTER
                for run in paragraph.runs:
                    set_run_font(run, size=font_size)


def build_document() -> None:
    doc = Document()
    section = doc.sections[0]
    section.page_width = Inches(8.5)
    section.page_height = Inches(11)
    section.top_margin = Inches(0.76)
    section.bottom_margin = Inches(0.72)
    section.left_margin = Inches(0.82)
    section.right_margin = Inches(0.82)

    styles = doc.styles
    normal = styles["Normal"]
    normal.font.name = "Times New Roman"
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    normal.font.size = Pt(10.5)
    normal.font.color.rgb = RGBColor(0, 0, 0)
    body = styles["Body Text"]
    body.font.name = "Times New Roman"
    body._element.rPr.rFonts.set(qn("w:eastAsia"), "宋体")
    body.font.size = Pt(10.5)
    body.font.color.rgb = RGBColor(0, 0, 0)

    for style_name, size in (("Title", 22), ("Heading 1", 15), ("Heading 2", 12), ("Heading 3", 10.5)):
        style = styles[style_name]
        style.font.name = "Arial"
        style._element.rPr.rFonts.set(qn("w:eastAsia"), "黑体")
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = RGBColor(0, 0, 0)

    title = doc.add_paragraph(style="Title")
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title.paragraph_format.space_before = Pt(82)
    title.paragraph_format.space_after = Pt(18)
    title_run = title.add_run("L 波段前向亮温计算模型公式和文献来源")
    set_run_font(title_run, latin="Arial", east_asia="黑体", size=22, bold=True)

    subtitle = doc.add_paragraph()
    subtitle.alignment = WD_ALIGN_PARAGRAPH.CENTER
    subtitle.paragraph_format.space_after = Pt(34)
    subtitle_run = subtitle.add_run("面向海冰厚度和冰上积雪联合反演的可审计实现说明")
    set_run_font(subtitle_run, latin="Arial", east_asia="黑体", size=12)
    subtitle_run.font.color.rgb = RGBColor(64, 64, 64)

    add_body(
        doc,
        "本文档逐步说明 tb_forward.m 中从海冰厚度、雪厚和冰型到 1.4 GHz 亮温的计算链，并把每个公式对应到论文、章节、公式号或表格。主复现对象是 2017 年 Data Synergy 论文；改进辐射模型论文、周璐博士论文和 Maaß 2013 用于补足其正文省略的物理细节。",
    )
    add_body(
        doc,
        "主要结论  当前默认配置应采用 20 层海冰和 1 层干雪。单层模型只用于敏感性基线。介电关系、层间辐射递推和顶层亮温分别来自不同文献层级；其中连续温盐廓线的离散递推属于本项目根据公开模型补全的实现，不等同于原作者源代码。",
        bold_lead="主要结论  ",
    )
    date_p = doc.add_paragraph()
    date_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    date_p.paragraph_format.space_before = Pt(46)
    r = date_p.add_run("整理日期  2026 年 10 月 18 日")
    set_run_font(r, size=9.5)
    r.font.color.rgb = RGBColor(89, 89, 89)
    doc.add_page_break()

    add_heading(doc, "一 模型范围和文献层级", 1)
    add_body(doc, "前向模型计算雪覆盖海冰在 L 波段的亮温。垂直介质依次为空气、干雪、多层海冰和半无限海水。2017 年联合反演论文第 2.1 节只把该过程概括为非线性函数，并说明采用改进多层模型；20 层海冰、温盐廓线和介电公式需要沿其引用链补齐。")
    add_formula(
        doc,
        "式 1  2017 联合反演论文使用的前向模型接口",
        "model_input",
        5.3,
        "Data Synergy 论文，第 2.1 节，式 (1)，PDF 第 4 页。该式定义输入输出关系，并未展开辐射传输细节。",
    )

    hierarchy_rows = [
        ("目标论文", "Data Synergy 2017", "第 2.1 节式 (1)；固定参数；表 1 八个理想场景", "决定复现目标和校验值"),
        ("改进模型", "Improving L-band radiation model 2017", "第 3 节式 (1) 至式 (3)；图 1；第 4.1 节", "决定 20 层冰、温盐廓线和多层方案"),
        ("详细汇总", "周璐博士论文", "第 3.2.2 节式 (3-1) 至式 (3-13)", "集中给出温度、盐度、介电与 TOA 亮温"),
        ("单层基线", "Maaß 2013", "第 2.1 节式 (1) 至式 (13)", "给出雪冰水体系、宇宙背景和稳态导热"),
        ("基础关系", "更早实验和模型文献", "Burke 1979；Kaleschke 2010；Cox and Weeks；Tiuri；Vant；Klein and Swift", "提供分层辐射和材料介电参数"),
    ]
    table = doc.add_table(rows=1 + len(hierarchy_rows), cols=4)
    populate_table(table, ["层级", "文献", "章节或公式", "在当前代码中的作用"], hierarchy_rows, [0.8, 1.7, 2.25, 2.1], 8.3)

    add_heading(doc, "二 计算流程", 1)
    for item in (
        "读取海冰厚度 h_i、雪厚 h_s、冰型、表面温度、海水温度与盐度。",
        "按层中点生成 FYI 或 MYI 盐度廓线，并用稳态导热计算雪层和冰层温度。",
        "由温度和盐度计算卤水体积分数，再求雪、冰和海水的复介电常数。",
        "对每个入射角和极化，计算 Fresnel 反射、层内衰减和逐层向上辐射。",
        "加入宇宙背景反射项和显式大气项，最后对 0 至 40 度及 H V 极化求平均。",
    ):
        add_bullet(doc, item)

    add_heading(doc, "三 参数值和来源", 1)
    parameter_rows = [
        ("海冰层数", "20", "Improving 论文第 3 节；博士论文第 3.2.2 节", "论文比较 1 20 100 层后采用 20"),
        ("积雪层数", "1", "博士论文第 3.2.2 节", "均匀干雪层"),
        ("频率", "1.4 GHz", "Data Synergy 第 2.1 节；Improving 第 2 节", "SMOS L 波段"),
        ("入射角", "0 到 40 度", "Improving 第 2 节；Zhou 2018 第 2.3 节", "1 度步长为代码离散选择"),
        ("极化", "H 和 V", "Improving 第 2 节；Zhou 2018 第 2.3 节", "两种极化算术平均"),
        ("表面温度", "负 30 摄氏度", "Data Synergy 第 2.1 节", "表 1 理想场景默认值"),
        ("海水温度", "负 1.8 摄氏度", "Data Synergy 第 2.1 节；Maaß 第 2.1 节", "冰底边界温度"),
        ("海水盐度", "33 g kg-1", "Data Synergy 第 2.1 节；Maaß 第 2.1 节", "Klein Swift 海水介电模型输入"),
        ("雪密度", "320 kg m-3", "Data Synergy 第 2.1 节；Improving 第 3 节", "前向函数换算为相对水密度 0.32"),
        ("雪热导率", "0.31 W m-1 K-1", "Improving 式 (2)；Maaß 第 2.1 节", "均匀雪层"),
        ("FYI 最小盐度", "1.5 ppt", "Improving 第 3 节；博士论文式 (3-5) 后", "代码采用下限保护"),
        ("宇宙背景", "2.7 K", "Maaß 第 2.1 节；Zhou 2018 附录 A3", "作为未发射部分的下行背景"),
        ("大气项", "0 K", "代码默认值", "论文只给负 0.36 到正 5.67 K 范围，未给逐场景值"),
        ("经验偏移", "0 K", "代码默认值", "不把验证偏差硬编码为物理项"),
        ("卤水适用上限", "70 permille", "Maaß 第 2.1 节；Kaleschke 第 2.2.1 节", "Vant 经验关系的标称有效域"),
        ("海冰密集度", "1", "Data Synergy 第 2.1 节的样本限制", "排除开阔水 水道和混合冰型"),
    ]
    table = doc.add_table(rows=1 + len(parameter_rows), cols=4)
    populate_table(table, ["参数", "当前值", "来源", "说明"], parameter_rows, [1.25, 1.25, 2.25, 2.1], 8.1)
    add_body(doc, "密度参数的使用需要区分。main_inversion.m 中的海水密度 1024 kg m-3 和整体海冰密度 915 kg m-3用于静力平衡反演，不直接进入 tb_forward.m。前向亮温函数使用雪密度 320 kg m-3，并在卤水体积计算中使用随温度变化的纯冰密度。")

    add_heading(doc, "四 温度廓线", 1)
    add_heading(doc, "四一 冰雪热导率", 2)
    add_formula(
        doc,
        "式 2  海冰和积雪热导率",
        "conductivity",
        5.4,
        "Improving 论文第 3 节式 (1) 和式 (2)；博士论文第 3.2.2 节式 (3-1) 和式 (3-2)；Maaß 第 2.1 节式 (9)。海冰式追溯至 Untersteiner 1964，采用参数也见 Yu and Rothrock 1996。",
    )
    add_body(doc, "代码用整根雪冰柱的平均温度代替尚未知的体平均冰温，这是 Maaß 2013 在式 (9) 后明确采用的简化。温度以摄氏度代入分母，等价于原式中绝对温度减 273.15。")

    add_heading(doc, "四二 稳态热阻和层温", 2)
    add_formula(
        doc,
        "式 3  稳态导热和线性层温",
        "heat_profile",
        4.6,
        "Maaß 第 2.1 节式 (4) 至式 (13)；Improving 第 3 节；博士论文第 3.2.2 节式 (3-3) 和式 (3-4)。",
    )
    add_body(doc, "代码先由串联热阻求热通量，再计算雪冰界面温度，并在每个数值层的中点取温度。雪层和冰层内部热导率各自视为常数，界面温度和热通量连续。")

    add_heading(doc, "五 盐度廓线", 1)
    add_heading(doc, "五一 一年冰盐度", 2)
    add_formula(
        doc,
        "式 4  FYI 厚度相关体盐度",
        "fyi_salinity",
        5.5,
        "博士论文第 3.2.2 节式 (3-5)；Improving 第 3 节说明采用厚度相关均匀盐度。原始观测关系引用 Cox and Weeks 1974。",
    )
    add_body(doc, "FYI 的所有冰层使用相同盐度。当前代码对经验式统一施加 1.5 ppt 下限；这是依据改进模型对厚冰盐度下限的处理。")

    add_heading(doc, "五二 多年冰盐度", 2)
    add_formula(
        doc,
        "式 5  MYI 层中点盐度廓线",
        "myi_salinity",
        5.4,
        "Improving 第 3 节式 (3) 和图 1；博士论文第 3.2.2 节式 (3-6) 和图 3.2。廓线形式依据 Schwarzacher 1959 的多年冰观测。",
    )
    add_formula(doc, "式 6  MYI 拟合参数", "myi_constants", 3.3, "Improving 第 3 节式 (3) 后给出的拟合值；博士论文式 (3-6) 后同样列出。")
    add_body(doc, "标准化深度从冰面 0 增加到冰底 1。该廓线反映 MYI 表层经历融化冲刷后盐度低、底部盐度较高的结构。")

    add_heading(doc, "六 卤水体积分数和介电常数", 1)
    add_heading(doc, "六一 卤水体积分数", 2)
    add_formula(
        doc,
        "式 7  温度相关纯冰密度和卤水体积分数",
        "brine",
        5.7,
        "Cox and Weeks 1983 式 (5) 至式 (7)、式 (14) 至式 (15)、表 I 和表 II；Kaleschke 2010 第 2.2.1 节式 (8) 与表 2；博士论文第 3.2.2 节式 (3-10) 至式 (3-13)。",
    )
    add_body(doc, "F1 和 F2 是温度的分段三次多项式。负 2 摄氏度以上采用 Kaleschke 2010 表 2 所列 Leppäranta Manninen 系数；负 2 至负 30 摄氏度采用 Cox and Weeks 1983 表 II，并在负 22.9 摄氏度处分段。代码最终将体积分数乘 1000，得到 permille。")

    add_heading(doc, "六二 干雪介电常数", 2)
    add_formula(
        doc,
        "式 8  干雪复介电常数",
        "snow_eps",
        6.1,
        "博士论文第 3.2.2 节式 (3-8)；Zhou 2018 附录 A3；原始实验拟合为 Tiuri et al. 1984。相对密度是雪密度除以 1000 kg m-3。",
    )
    add_body(doc, "当前模型只处理干雪，不包含湿度、风砌雪、深霜或垂直各向异性。该限制来自 Improving 论文第 3 节和博士论文第 3.2.2 节。")

    add_heading(doc, "六三 海冰介电常数", 2)
    add_formula(
        doc,
        "式 9  海冰复介电常数",
        "ice_eps",
        3.7,
        "Vant et al. 1978；Kaleschke 2010 第 2.2.1 节式 (7) 和表 1；博士论文第 3.2.2 节式 (3-9)。",
    )
    coeff_rows = [
        ("FYI", "3.10", "0.0084", "0.037", "0.00445"),
        ("MYI", "3.10", "0.0084", "0.003", "0.00435"),
    ]
    table = doc.add_table(rows=1 + len(coeff_rows), cols=5)
    populate_table(table, ["冰型", "a1", "a2", "a3", "a4"], coeff_rows, [1.4, 1.2, 1.2, 1.2, 1.2], 8.8)
    add_body(doc, "这些系数对应 1.4 GHz。Vant 关系标称适用于卤水体积分数小于 70 permille。代码保留原始最大值用于诊断，同时把进入介电关系的数值限制到 70 permille。")

    add_heading(doc, "六四 海水介电常数", 2)
    add_formula(
        doc,
        "式 10  Klein Swift 单 Debye 海水模型",
        "water_eps",
        6.0,
        "Klein and Swift 1977；Maaß 第 2.1 节；博士论文第 3.2.2 节；Zhou 2018 附录 A3。当前代码展开了原文给出的静态介电常数、弛豫时间和电导率温盐多项式。",
    )
    add_body(doc, "代码采用无穷频率介电常数 4.9，并把海水温度和盐度代入 Klein Swift 多项式。复数符号最终转换为与其余模型一致的正损耗虚部约定。")

    add_heading(doc, "七 Fresnel 边界和分层辐射", 1)
    add_heading(doc, "七一 边界反射", 2)
    add_formula(
        doc,
        "式 11  H 和 V 极化 Fresnel 功率反射率",
        "fresnel",
        6.25,
        "Kaleschke 2010 第 2.1 至第 2.2 节；Maaß 第 2.1 节；博士论文第 3.2.2 节引用 Ulaby 等微波遥感教材。",
    )
    add_body(doc, "空气介电常数取 1。折射角通过横向波数守恒隐含在 q 中，因此公式可以直接使用外部入射角。")

    add_heading(doc, "七二 层内衰减", 2)
    add_formula(
        doc,
        "式 12  单程功率透射率",
        "attenuation",
        5.75,
        "Burke et al. 1979 的分层非相干辐射框架；Maaß 第 2.1 节。Kaleschke 2010 第 2.2 节式 (3) 至式 (4) 给出相应传播常数和衰减系数。",
    )

    add_heading(doc, "七三 向上辐射递推", 2)
    add_formula(
        doc,
        "式 13  当前代码使用的一阶非相干层递推",
        "recursion",
        4.8,
        "物理框架来自 Burke et al. 1979、Maaß 第 2.1 节和 Zhou 2018 附录 A1 与 A3。A U E 三个递推式是本项目为可运行实现写出的离散形式，公开论文没有提供完全相同的逐行公式或源代码。",
    )
    add_body(doc, "递推从半无限海水向上进行。U 表示热辐射亮温，E 表示单位底边界输入所对应的有效发射率权重。雪冰界面和冰水界面计算 Fresnel 反射；同一种材料内部的相邻数值层被视为连续廓线的积分单元，不额外施加人造界面反射。")
    add_formula(doc, "式 14  顶界面输出", "top_surface", 4.5, "Burke 分层辐射框架和 Maaß 第 2.1 节。顶界面为空气与最上层介质的 Fresnel 边界。")

    add_heading(doc, "八 顶层亮温和 SMOS 平均", 1)
    add_heading(doc, "八一 宇宙背景和大气项", 2)
    add_formula(
        doc,
        "式 15  纯冰覆盖条件下的顶层亮温",
        "toa",
        6.0,
        "Maaß 第 2.1 节的海冰和开阔水亮温式；博士论文第 3.2.2 节式 (3-7)；Zhou 2018 附录 A3 式 (A1)。这里已按海冰密集度 c=1 化简。",
    )
    add_body(doc, "宇宙背景固定为 2.7 K。Zhou 2018 给出大气贡献约为负 0.36 至正 5.67 K，但没有提供每个理想场景的确定值，因此代码默认大气项为 0。经验偏移同样默认 0；论文中的拟合截距或验证偏差不能直接解释为纯大气校正。")

    add_heading(doc, "八二 角度和极化平均", 2)
    add_formula(
        doc,
        "式 16  SMOS 日产品约定下的平均亮温",
        "angular_mean",
        5.3,
        "Improving 第 2 节说明产品是 0 至 40 度、H V 极化平均；Zhou 2018 第 2.3 节明确按同一协议模拟。代码使用 1 度等间距采样，这是数值实现选择。",
    )

    add_heading(doc, "九 代码实现与论文原式的对应", 1)
    mapping_rows = [
        ("apply_defaults", "层数 频率 角度和固定参数", "Improving 第 2 至第 3 节；Data Synergy 第 2.1 节；Zhou 2018 附录 A"),
        ("vertical_profiles", "FYI MYI 盐度及线性温度", "博士论文式 (3-1) 至式 (3-6)；Improving 式 (1) 至式 (3)"),
        ("brine_volume_permille", "纯冰密度及 F1 F2 分段卤水体积", "Cox and Weeks 1983；Kaleschke 式 (8) 和表 2"),
        ("snow_permittivity", "干雪介电常数", "博士论文式 (3-8)；Tiuri 1984"),
        ("ice_permittivity", "Vant 海冰介电关系和冰型系数", "博士论文式 (3-9)；Kaleschke 式 (7) 和表 1"),
        ("seawater_permittivity", "Klein Swift 海水 Debye 模型", "Klein and Swift 1977；Maaß 第 2.1 节"),
        ("fresnel_power", "H V 极化边界反射", "Kaleschke 第 2.1 至第 2.2 节；Ulaby 教材"),
        ("burke_stack", "层衰减 自发辐射和向上递推", "Burke 1979；Maaß 第 2.1 节；本项目离散实现"),
        ("forward_scalar", "宇宙背景 大气项和多角度平均", "Zhou 2018 式 (A1)；Improving 第 2 节"),
    ]
    table = doc.add_table(rows=1 + len(mapping_rows), cols=3)
    populate_table(table, ["代码位置", "计算内容", "对应文献"], mapping_rows, [1.6, 2.35, 2.85], 8.25)

    add_heading(doc, "十 实现选择和适用边界", 1)
    add_bullet(doc, "连续温盐廓线被离散为数值层。同一种材料内部不增加 Fresnel 界面；否则微小介电差的反射会随层数累积，无法得到合理的层数收敛。")
    add_bullet(doc, "Burke 1979 和 Maaß 2013 给出分层非相干模型的物理框架，但没有公开与当前函数完全相同的源代码。因此 burke_stack 属于可审计重建。")
    add_bullet(doc, "当前模型假设干雪、平滑界面、冬季状态、单一冰型和海冰密集度为 1，不模拟湿雪、融池、表面粗糙度、开阔水或冰间水道面积混合。")
    add_bullet(doc, "Vant 海冰介电关系在卤水体积分数低于 70 permille 时有效。代码对介电计算输入限幅，并在 diagnostics 中报告原始最大值和越界标志。")
    add_bullet(doc, "大气亮温和经验偏移保持显式参数且默认 0。若以后用观测标定，应单独保存标定依据，不能把拟合偏差直接写成物理常数。")

    add_heading(doc, "十一 论文场景回归检验", 1)
    add_body(doc, "validate_tb_forward.m 的第一组回归检验逐项采用 Data Synergy 论文表 1 的八个理想场景。下表给出论文亮温、当前模型亮温和差值。")
    validation_rows = [
        ("I", "FYI", "0.50", "0.05", "231.5", "230.726", "-0.774"),
        ("II", "FYI", "1.00", "0.03", "238.2", "237.796", "-0.404"),
        ("III", "FYI", "1.50", "0.10", "243.7", "242.582", "-1.118"),
        ("IV", "FYI", "2.50", "0.25", "246.5", "244.782", "-1.718"),
        ("V", "MYI", "1.50", "0.15", "243.8", "242.374", "-1.426"),
        ("VI", "MYI", "2.50", "0.15", "248.9", "246.415", "-2.485"),
        ("VII", "MYI", "3.00", "0.35", "252.0", "248.964", "-3.036"),
        ("VIII", "MYI", "5.00", "0.40", "250.4", "246.859", "-3.541"),
    ]
    table = doc.add_table(rows=1 + len(validation_rows), cols=7)
    populate_table(table, ["场景", "冰型", "h i m", "h s m", "论文 TB K", "模型 TB K", "差值 K"], validation_rows, [0.62, 0.65, 0.68, 0.68, 1.05, 1.05, 0.85], 7.9)
    add_body(doc, "八个场景的 RMSE 为 2.0895 K，平均偏差为负 1.8128 K。另一个测试在负 35 摄氏度、1 m 冰和 0.1 m 雪条件下得到 FYI 与 MYI 亮温差 11.28 K，接近改进模型论文报告的约 12 K。5 m 厚 MYI 的 20 层与 100 层差小于 0.1 K；薄冰仍有更强的层数敏感性。")

    add_heading(doc, "十二 参考文献和链接", 1)
    add_body(doc, "以下前四项为用户已提供文献，按要求不重复附链接。其余文献提供 DOI 或出版机构页面。")

    provided_refs = [
        "Zhou L 等  Data Synergy between Altimetry and L-Band Passive Microwave Remote Sensing for the Retrieval of Sea Ice Parameters  A Theoretical Study of Methodology  Remote Sensing 2017",
        "Zhou L 等  Improving L-band radiation model and representation of small-scale variability to simulate brightness temperature of sea ice  International Journal of Remote Sensing 2017",
        "周璐  北极海冰厚度和冰上积雪的联合反演及相关研究  博士论文",
        "Maaß N 等  Snow thickness retrieval over thick Arctic sea ice using SMOS satellite data  The Cryosphere 2013",
    ]
    for idx, ref in enumerate(provided_refs, start=1):
        p = doc.add_paragraph(style="List Number")
        p.paragraph_format.space_after = Pt(3)
        set_run_font(p.add_run(ref), size=9.3)

    external_refs = [
        ("Kaleschke L 等  A sea-ice thickness retrieval model for 1.4 GHz radiometry and application to airborne measurements over low salinity sea-ice  The Cryosphere 2010", "https://doi.org/10.5194/tc-4-583-2010"),
        ("Zhou L 等  On the retrieval of sea ice thickness and snow depth using concurrent laser altimetry and L-band remote sensing data  The Cryosphere 2018", "https://doi.org/10.5194/tc-12-993-2018"),
        ("Burke W J  Schmugge T  Paris J F  Comparison of 2.8- and 21-cm microwave radiometer observations over soils with emission model calculations  1979", "https://doi.org/10.1029/JC084iC01p00287"),
        ("Cox G F N  Weeks W F  Salinity Variations in Sea Ice  1974", "https://doi.org/10.3189/S0022143000023418"),
        ("Cox G F N  Weeks W F  Equations for Determining the Gas and Brine Volumes in Sea-Ice Samples  1983", "https://doi.org/10.1017/S0022143000008364"),
        ("Schwarzacher W  Pack Ice Studies in the Arctic Ocean  1959", "https://doi.org/10.1029/JZ064i012p02357"),
        ("Tiuri M 等  The complex dielectric constant of snow at microwave frequencies  1984", "https://doi.org/10.1109/JOE.1984.1145645"),
        ("Vant M R  Ramseier R O  Makios V  The complex dielectric constant of sea ice at frequencies in the range 0.1 to 40 GHz  1978", "https://doi.org/10.1063/1.325018"),
        ("Klein L A  Swift C T  An improved model for the dielectric constant of sea water at microwave frequencies  1977", "https://doi.org/10.1109/TAP.1977.1141539"),
        ("Untersteiner N  Calculations of temperature regime and heat budget of sea ice in the central Arctic  1964", "https://doi.org/10.1029/JZ069i022p04755"),
        ("Yu Y  Rothrock D A  Thin ice thickness from satellite thermal imagery  1996", "https://doi.org/10.1029/96JC02242"),
        ("Wadhams P 等  Relationship between sea ice freeboard and draft in the Arctic Basin and implications for ice thickness monitoring  1992", "https://doi.org/10.1029/92JC02014"),
        ("Warren S G 等  Snow Depth on Arctic Sea Ice  Journal of Climate 1999", "https://doi.org/10.1175/1520-0442(1999)012%3C1814:SDOASI%3E2.0.CO;2"),
    ]
    start = len(provided_refs) + 1
    for idx, (ref, url) in enumerate(external_refs, start=start):
        p = doc.add_paragraph(style="List Number")
        p.paragraph_format.space_after = Pt(3)
        set_run_font(p.add_run(ref + "  "), size=9.3)
        add_hyperlink(p, "在线链接", url)

    add_heading(doc, "十三 使用建议", 1)
    add_body(doc, "论文复现应以 20 层海冰和 1 层雪作为默认方案。单层结果只用于确认模型演化和灵敏度。启动全量反演前，应先决定是否接受当前表 1 约 2.09 K 的误差，并继续检查薄冰层数收敛和高卤水体积分数场景。")

    doc.core_properties.title = "L 波段前向亮温计算模型公式和文献来源"
    doc.core_properties.subject = "海冰厚度和冰上积雪联合反演前向模型"
    doc.core_properties.author = "Codex"
    doc.core_properties.keywords = "SMOS L-band sea ice snow brightness temperature"
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    doc.save(OUTPUT)


if __name__ == "__main__":
    build_document()
    print(OUTPUT)
