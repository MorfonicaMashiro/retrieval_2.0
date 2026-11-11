"""从同一组图元生成可编辑 draw.io 与预览 PNG。"""
from pathlib import Path
import xml.etree.ElementTree as ET
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parent
STEM = ROOT / "联合反演实现与实验分支"
W, H = 2460, 3670
FONT_PATH = r"C:\Windows\Fonts\msyh.ttc"
COLORS = {
    "navy": "#173B63", "blue": "#EAF3FC", "blue_line": "#5B8DB8",
    "teal": "#E7F5F1", "teal_line": "#469587",
    "orange": "#FFF1E1", "orange_line": "#C58A48",
    "purple": "#F2EDFA", "purple_line": "#8D75AD",
    "gray": "#F2F4F7", "gray_line": "#84909C",
    "white": "#FFFFFF", "ink": "#243341", "muted": "#526577",
}

items = []
edges = []

def normalize_display(value):
    """保留公式含义，同时避开常见中文字库缺失的数学下标字形。"""
    replacements = {
        "kg·m⁻³": "kg/m^3", "W·m⁻¹·K⁻¹": "W/(m·K)",
        "10⁻⁵": "1e-5", "sin²": "sin^2", "R²": "R^2",
        "d⁻ᵖ": "d^(-p)", "q̄": "q_mean", "⌊": "floor(", "⌋": ")",
        "Uₙ₊₁": "U_(n+1)", "Uⱼ₊₁": "U_(j+1)",
        "₀": "_0", "₁": "_1", "₂": "_2", "₃": "_3", "ₙ": "_n",
        "ⱼ": "_j", "ₖ": "_k", "ₘ": "_m", "ₚ": "_p", "ᵢ": "_i",
        "⁻": "-", "¹": "^1", "²": "^2", "³": "^3",
    }
    for before, after in replacements.items():
        value = value.replace(before, after)
    return value

def box(key, x, y, w, h, title, lines, fill, stroke, fs=23, title_fs=25):
    items.append(dict(key=key, x=x, y=y, w=w, h=h,
                      title=normalize_display(title),
                      lines=[normalize_display(line) for line in lines],
                      fill=fill, stroke=stroke, fs=fs, title_fs=title_fs))

def edge(key, coords, color="#63809A", dashed=False, arrow=True):
    edges.append(dict(key=key, coords=coords, color=color, dashed=dashed, arrow=arrow))

def section(x, y, label, w=1400):
    box("section_" + label, x, y, w, 62, label, [], COLORS["navy"], COLORS["navy"], 23, 27)

section(80, 145, "① 同日数据获取与 EASE 配准")
box("smos_raw", 80, 232, 680, 216, "SMOS 逐日亮温", [
    "ICDC L3B v3，2012/2015 OIB 飞行同日", "原生 PSP/EPSG:3413，896×608；TB、nPair",
    "读取 NetCDF；−999 → NaN；不做季节平均"], COLORS["blue"], COLORS["blue_line"])
box("oib_raw", 800, 232, 680, 216, "OIB 逐航次观测", [
    "2012 IDCSI4 L4；2015 NSIDC-0708 Quick Look", "lat/lon、mean_fb、snow_depth、thickness、冰型",
    "−99999 → NaN；同日文件配对"], COLORS["blue"], COLORS["blue_line"])
box("smos_ease", 80, 495, 680, 230, "SMOS → 原版北极 EASE", [
    "目标 1441×1441；Δ=12 533.76 m", "目标中心 (x,y) → 经纬度 → 源 PSP",
    "最近源像元复制 TB；nPair≤0 保持缺测",
    "逐中心变换，不能仅把整图旋转 45°"], COLORS["teal"], COLORS["teal_line"])
box("oib_ease", 800, 495, 680, 230, "OIB 点 → EASE 行列", [
    "经纬度直接投影；数值不插值", "col=⌊(x−xMin)/Δ⌋+1",
    "row=⌊(yMax−y)/Δ⌋+1", "FBice,p=mean_fb,p−snow_depth,p"], COLORS["teal"], COLORS["teal_line"])
box("collocate", 80, 775, 1400, 260, "同日 3×3 配准与共同样本集合", [
    "每个有效 SMOS 中心格取本格及周围 8 个 EASE 格内的 OIB 原始点（约 37.5×37.5 km）",
    "四量有限；开水/薄冰/灰冰比例均为 0；FYI/MYI 单一冰型；nOIB>0；TB 有效",
    "对同一批 M 个点：q̄=(1/M)Σqₚ，q∈{FBsnow, hs, FBice, hi}；不先求九格均值",
    "中心格 stride=1，可重叠；不额外排除负 FBice（求解器再检查）"], COLORS["blue"], COLORS["blue_line"], 22)
box("input", 80, 1082, 1400, 165, "一条配对的反演输入与验证参照", [
    "输入：TBobs、冰型、FBice 或 FBsnow；固定 ρw/ρi/ρs=1024/915/320 kg·m⁻³",
    "OIB 平均 hi、hs 只作验证；FBice 是雪面干舷减雪深的代理，并非独立雷达观测",
    "同一配对独立求解；记录 day/row/col/nOIB、无根与双根状态"], COLORS["orange"], COLORS["orange_line"], 22)

section(80, 1285, "② F(hi, hs, 冰型)：1.4 GHz 多层前向亮温")
box("profile", 80, 1372, 1400, 245, "温度和盐度廓线（20 冰层 + 1 雪层）", [
    "稳态热阻 R=hs/ks+hi/ki；q=(Tw−Tsurf)/R；T界=Tsurf+q·hs/ks",
    "雪/冰各层取中点 T；Tsurf=−30 °C，Tw=−1.8 °C，ks=0.31 W·m⁻¹·K⁻¹",
    "ki=2.034+0.13·mean(S)/Tbulk,C；Tbulk,C≈(Tsurf+Tw)/2",
    "FYI: S=max[1.5, 6.08e^(−5.81hi)+7.409e^(−0.5228hi)] ppt，各层相同",
    "MYI: zₖ=(k−½)/N；Sₖ=1.6[1−cos(πzₖ^(0.407/(zₖ+0.573)))] ppt"], COLORS["teal"], COLORS["teal_line"], 22)
box("dielectric", 80, 1638, 1400, 252, "卤水体积与雪—冰—海水复介电常数", [
    "Cox–Weeks: Vb=1000·ρ纯冰(T)·S/F₁(T) (‰)；F₁ 随温度分段；>70‰ 诊断标记",
    "Vant: εice=a₁+a₂Vb+i(a₃+a₄Vb)；FYI/MYI 使用各自 1.4 GHz 系数",
    "Tiuri 干雪: ε′=1+0.7ρr+0.7ρr²；ε″=C(f,ρr)·exp(0.036T)",
    "Klein–Swift 海水: 先算 Debye εraw，再取 Re(εraw)+i|Im(εraw)|；Sw=33 ppt",
    "默认 Vb 不截断、界面用相邻层中点 ε；层源项始终取层中点 T"], COLORS["teal"], COLORS["teal_line"], 22)
box("optics", 80, 1928, 1400, 197, "Fresnel 反射、层透射与有效发射率", [
    "qⱼ=√(εⱼ−sin²θ₀)；rH=(q₁−q₂)/(q₁+q₂)；rV=(ε₂q₁−ε₁q₂)/(ε₂q₁+ε₁q₂)",
    "R=|r|²；τⱼ=exp[−2k₀|Im(qⱼ)|dⱼ]，k₀=2πf/c；同材数值层界面 R=0",
    "以单位温度源重复同一递推得 e_eff；背景反射项为 (1−e_eff)×2.7 K"], COLORS["teal"], COLORS["teal_line"], 22)
box("rt", 80, 2163, 1400, 240, "Burke 一阶非相干亮温递推：海水 → 冰 → 雪 → 空气", [
    "底边界 Uₙ₊₁=Tw(K)；自下而上对层 j：",
    "Aⱼ=(1−τⱼ)(1+Rⱼτⱼ)",
    "Uⱼ=Tⱼ(K)Aⱼ+(1−Rⱼ)τⱼUⱼ₊₁",
    "顶层：TBθ,p=(1−Rtop)U₁+(1−e_eff)·2.7 K；默认大气项/经验偏移均为 0"], COLORS["teal"], COLORS["teal_line"], 22)
box("average", 80, 2439, 1400, 137, "模拟亮温的观测口径", [
    "F(hi,hs,type)=mean_{θ=0:1:40°, p∈{H,V}} TBθ,p；先各角度/极化计算再算术平均",
    "纯冰覆盖样本：不加开阔水/水道面积混合；SMOS TB 与 F 在配对尺度比较"], COLORS["orange"], COLORS["orange_line"], 22)

section(80, 2612, "③ 两种干舷路线：扫描 + 二分求根")
box("ice_route", 80, 2699, 680, 279, "路线 A｜TB + 冰面干舷 FBice", [
    "候选 hs∈[0,1 m]；4 cm 初扫",
    "hi=(ρw·FBice+ρs·hs)/(ρw−ρi)",
    "计算 g(hs)=F(hi,hs,type)−TBobs",
    "保留满足 0≤hi≤10 m、FBice≥0 的点",
    "扫描命中或变号区间交给二分"], COLORS["orange"], COLORS["orange_line"], 22)
box("snow_route", 800, 2699, 680, 279, "路线 B｜TB + 雪面干舷 FBsnow", [
    "候选 hs∈[0,min(1 m, FBsnow)]；1 cm 全扫",
    "hi=[ρw·FBsnow−(ρw−ρs)hs]/(ρw−ρi)",
    "FBice=FBsnow−hs≥0；0≤hi≤10 m",
    "每个变号区间都二分；保留全部根",
    "验证默认取较小 hs 根，第二支另存"], COLORS["orange"], COLORS["orange_line"], 22)
box("solver", 80, 3020, 1400, 180, "数值终止与解记录", [
    "对 g(hs) 异号区间取中点并更新边界；扫描点 |g|≤0.01 K 也可直接命中",
    "停止条件：|g|≤0.01 K 或 hs 区间宽≤10⁻⁵ m；最多 60 次；最终再次校验 TB 残差",
    "存 hi、hs、TBsim、残差、根数/无根状态；多根不强行合并为单解"], COLORS["blue"], COLORS["blue_line"], 22)
box("metrics", 80, 3238, 1400, 270, "④ 误差评定与精度保留", [
    "与同一配对 OIB 平均 hi/hs 比较；TB 验证值另算 F(hiOIB,hsOIB,type)",
    "有效成对样本上：bias=mean(反演−参照)，RMSE=√mean[(反演−参照)²]",
    "自由斜率最小二乘 y=ax+b、相关系数 r 与 R²=r²；保留 n、无解率、双解率",
    "同时报告全部/逐年/前 5% nOIB/|TBOIB−TBobs|≤1.5 K 诊断子集",
    "用原始双精度与具名缓存；不按打印位数提前舍入；1.5 K 子集不能作真实反演筛选器"], COLORS["blue"], COLORS["blue_line"], 22)

section(1560, 145, "实验分支与接入位置", 820)
box("tests_data", 1560, 232, 820, 435, "A 数据与空间算子｜已执行/已实现", [
    "基线：2012/2015，stride=1，3×3 重叠；6,542 配对",
    "扩展：2013/2014 同日配准；四年共 16,364 配对",
    "stride=3、phase=[2,2]：3×3 不重叠，733 配对",
    "参数先均值 A：F(mean hi, mean hs)",
    "逐点 TB 均值 B：mean F(hiₘ,hsₘ)",
    "加权参数 Aw：F(Σw hi/Σw,Σw hs/Σw)",
    "加权逐点 Bw：Σw F(hiₘ,hsₘ)/Σw",
    "权重：中心反距离 w=d⁻ᵖ（正式 p=1）；",
    "或 2.5 km 等面积子网格；逐点反演已实现"], COLORS["purple"], COLORS["purple_line"], 20, 24)
box("tests_physics", 1560, 708, 820, 455, "B 前向物理与离散｜已执行/可配置", [
    "表 1 八场景 TB 对照；FYI/MYI、温盐/界面审查",
    "Cox–Weeks F₁-only 与旧 F₂ 闭合口径审计",
    "Vb>70‰：不截断基线 vs 70‰ 截断敏感性",
    "冰层 N=1/20/100；20→100 层差单列报告",
    "物理界面 ε：midpoint vs boundary",
    "辐射递推：Burke vs upperReflection",
    "宇宙背景/大气/偏移开关仅作具名敏感性",
    "雪层数可配置；目标固定层厚度试验尚未实现",
    "不同物理口径需独立 variantTag 与结果文件"], COLORS["purple"], COLORS["purple_line"], 20, 24)
box("tests_solver", 1560, 1204, 820, 405, "C 反演与闭环｜已执行/已实现", [
    "八理想场景 × 两条路线：模型自闭环",
    "论文表 1 → 表 2：独立检验可解性/双解",
    "扫描步长 4 cm / 1 cm、容差、无根状态",
    "逐点干舷求根：每个候选 hs 对所有点重算 TB",
    "雪面逐点路线：局部 hs=min(候选,FBsnowₘ)",
    "淹没权重占比>0.5 则候选无效",
    "等权算子 2 与等面积加权算子 4 均有入口",
    "逐点路线的 hi/hs 参照使用同一组空间权重"], COLORS["purple"], COLORS["purple_line"], 20, 24)
box("tests_validation", 1560, 1650, 820, 376, "D 评定分支｜已执行", [
    "TB 全样本、nOIB 前 5% 与逐年指标",
    "冰厚/雪深 × 两种干舷路线，各报 n/r/R²/bias/RMSE",
    "|TBOIB−SMOS|≤1.5 K 诊断子集另报",
    "FYI/MYI、算子 A/B 差值与中心化 RMSE 对照",
    "保存第二支雪面解；统计无根/双根及超 70‰ 比例",
    "基线输出冻结；扩展保持同批窗口和独立文件名",
    "OIB thickness 源于几何/静力平衡，非完全独立真值"], COLORS["purple"], COLORS["purple_line"], 20, 24)
box("tests_next", 1560, 2066, 820, 351, "E 后续可接入｜计划，勿标为已跑", [
    "雪层数 1/2/5/10 与目标层厚度收敛",
    "逐点 KT19_surf 温度、SMOS 质量/RFI/不确定度",
    "不同 nOIB 门槛、足迹响应与协变分析",
    "参数扰动、输入误差传播与样本独立性审查",
    "CryoSat-2、MODIS、海冰密集度、ERA5 等",
    "属于卫星尺度/周璐论文后续数据链，尚未并入",
    "这张 2017 OIB–SMOS 基线反演链"], COLORS["gray"], COLORS["gray_line"], 20, 24)
box("legend", 1560, 2461, 820, 173, "读图约定", [
    "蓝：观测/配准与评定；青：前向物理", "橙：反演输入/两路求解；紫：实验分支",
    "实线：主流程；虚线：测试从对应环节分出"], COLORS["gray"], COLORS["gray_line"], 20, 24)
box("note", 1560, 2677, 820, 331, "复现边界", [
    "本图对应当前 2017 Data Synergy MATLAB 复现。",
    "本地 SMOS 是 ICDC L3B PSP；论文写 EASE/L4，",
    "其重采样算法与完整原作者离散代码未公开。",
    "最近邻 EASE、Burke 递推及 Vb 外推均为",
    "本项目可审计的实现选择；不能称逐位复现。",
    "正式基线：Burke / midpoint / Vb 不截断 /",
    "20 冰层 + 1 雪层 / 2.7 K 宇宙背景。"], COLORS["gray"], COLORS["gray_line"], 20, 24)
box("source", 1560, 3045, 820, 463, "主要代码入口与结果口径", [
    "main.m → read_smos/read_oib → regrid_smos_to_ease",
    "→ PSP2EASE / collocate_day.m",
    "main_inversion.m → invert_pair.m",
    "→ tb_forward.m → burke_stack.m",
    "retrieval_config.m 统一科学参数",
    "validate_tb_forward.m、validate_inversion_closure.m",
    "analysis/tb_spatial_averaging/：四种算子与逐点反演",
    "基线结果：inversion_2012_2015_",
    "paper_baseline_v1_uncapped.*", "日期：2026-10-30"], COLORS["gray"], COLORS["gray_line"], 20, 24)

edge("e_smos", [(420,449),(420,493)])
edge("e_oib", [(1140,449),(1140,493)])
edge("e_smos_join", [(420,726),(420,750),(780,750),(780,774)])
edge("e_oib_join", [(1140,726),(1140,750),(780,750)], arrow=False)
for key, y1, y2 in [
    ("e_coll_input",1036,1081),("e_input_profile",1248,1284),
    ("e_section_profile",1348,1371),
    ("e_prof_diel",1618,1637),("e_diel_optics",1891,1927),
    ("e_optics_rt",2126,2162),("e_rt_avg",2404,2438),
    ("e_avg_routes",2577,2611)]:
    edge(key, [(780,y1),(780,y2)])
edge("e_branch_ice", [(780,2675),(780,2685),(420,2685),(420,2698)])
edge("e_branch_snow", [(780,2685),(1140,2685),(1140,2698)])
edge("e_ice_solve", [(420,2979),(420,3004),(780,3004),(780,3019)])
edge("e_snow_solve", [(1140,2979),(1140,3004),(780,3004)], arrow=False)
edge("e_solve_eval", [(780,3201),(780,3237)])
for key, y1, y2 in [
    ("t_data",875,450),("t_phys",1790,930),("t_solver",2890,1410),
    ("t_eval",3360,1830)]:
    edge(key, [(1481,y1),(1518,y1),(1518,y2),(1559,y2)], "#9B7BBD", True)

def drawio():
    mxfile = ET.Element("mxfile", {"host":"app.diagrams.net", "agent":"Codex", "version":"24.7.17", "pages":"1"})
    diagram = ET.SubElement(mxfile,"diagram",{"id":"seaice-flow","name":"联合反演总流程与实验分支"})
    model = ET.SubElement(diagram,"mxGraphModel",{"dx":str(W),"dy":str(H),"grid":"0", "guides":"1", "page":"1", "pageScale":"1", "pageWidth":str(W), "pageHeight":str(H), "math":"0", "shadow":"0"})
    root = ET.SubElement(model,"root")
    ET.SubElement(root,"mxCell",{"id":"0"}); ET.SubElement(root,"mxCell",{"id":"1","parent":"0"})
    title = ET.SubElement(root,"mxCell",{"id":"title","value":"海冰厚度—雪深联合反演：已实现流程与实验分支", "style":"text;html=1;strokeColor=none;fillColor=none;fontSize=38;fontStyle=1;fontColor=#173B63;fontFamily=Microsoft YaHei;align=center;verticalAlign=middle;", "vertex":"1","parent":"1"})
    ET.SubElement(title,"mxGeometry",{"x":"80","y":"35","width":"2300","height":"84","as":"geometry"})
    for item in items:
        is_section = item["key"].startswith("section_")
        if is_section:
            value = item["title"]
        else:
            value = "<b>" + item["title"] + "</b><br>" + "<br>".join(item["lines"])
        style = (f"rounded={0 if is_section else 1};arcSize=10;whiteSpace=wrap;html=1;"
                 f"fillColor={item['fill']};strokeColor={item['stroke']};strokeWidth={1 if is_section else 2};"
                 f"fontSize={item['fs']};fontColor={'#FFFFFF' if is_section else COLORS['ink']};"
                 "fontFamily=Microsoft YaHei;align=left;verticalAlign=middle;spacingLeft=18;spacingRight=12;spacingTop=10;")
        cell = ET.SubElement(root,"mxCell",{"id":item["key"],"value":value,"style":style,"vertex":"1","parent":"1"})
        ET.SubElement(cell,"mxGeometry",{"x":str(item["x"]),"y":str(item["y"]),"width":str(item["w"]),"height":str(item["h"]),"as":"geometry"})
    for ed in edges:
        style = ("edgeStyle=orthogonalEdgeStyle;rounded=0;html=1;"
                 f"strokeColor={ed['color']};strokeWidth=3;"
                 f"dashed={1 if ed['dashed'] else 0};"
                 f"endArrow={'block' if ed['arrow'] else 'none'};endFill=1;endSize=8;")
        cell = ET.SubElement(root,"mxCell",{"id":ed["key"],"value":"","style":style,"edge":"1","parent":"1"})
        geom = ET.SubElement(cell,"mxGeometry",{"relative":"1","as":"geometry"})
        ET.SubElement(geom,"mxPoint",{"x":str(ed["coords"][0][0]),"y":str(ed["coords"][0][1]),"as":"sourcePoint"})
        ET.SubElement(geom,"mxPoint",{"x":str(ed["coords"][-1][0]),"y":str(ed["coords"][-1][1]),"as":"targetPoint"})
        if len(ed["coords"]) > 2:
            pts = ET.SubElement(geom,"Array",{"as":"points"})
            for x,y in ed["coords"][1:-1]: ET.SubElement(pts,"mxPoint",{"x":str(x),"y":str(y)})
    ET.indent(mxfile, space="  ")
    ET.ElementTree(mxfile).write(STEM.with_suffix(".drawio"), encoding="utf-8", xml_declaration=True)

def rgba(hexvalue): return tuple(int(hexvalue[i:i+2],16) for i in (1,3,5))

def preview():
    img=Image.new("RGB",(W,H),"#FBFCFE")
    d=ImageDraw.Draw(img)
    ftitle=ImageFont.truetype(FONT_PATH,38)
    d.text((W/2,75),"海冰厚度—雪深联合反演：已实现流程与实验分支",font=ftitle,fill=COLORS["navy"],anchor="mm")
    for ed in edges:
        pts=ed["coords"]
        for a,b in zip(pts,pts[1:]):
            if ed["dashed"]:
                length=((b[0]-a[0])**2+(b[1]-a[1])**2)**0.5
                n=max(1,int(length/16))
                for i in range(0,n,2):
                    p=(a[0]+(b[0]-a[0])*i/n,a[1]+(b[1]-a[1])*i/n)
                    q=(a[0]+(b[0]-a[0])*min(i+1,n)/n,a[1]+(b[1]-a[1])*min(i+1,n)/n)
                    d.line([p,q],fill=ed["color"],width=3)
            else: d.line([a,b],fill=ed["color"],width=3)
        if ed["arrow"]:
            x,y=pts[-1]; a=pts[-2]
            if y>a[1]: poly=[(x,y),(x-8,y-13),(x+8,y-13)]
            elif x>a[0]: poly=[(x,y),(x-13,y-8),(x-13,y+8)]
            elif x<a[0]: poly=[(x,y),(x+13,y-8),(x+13,y+8)]
            else: poly=[(x,y),(x-8,y+13),(x+8,y+13)]
            d.polygon(poly,fill=ed["color"])
    for it in items:
        x,y,w,h=it["x"],it["y"],it["w"],it["h"]
        sec=it["key"].startswith("section_")
        d.rounded_rectangle((x,y,x+w,y+h),radius=0 if sec else 14,fill=it["fill"],outline=it["stroke"],width=2)
        f=ImageFont.truetype(FONT_PATH,it["title_fs"])
        if sec:
            d.text((x+18,y+h/2),it["title"],font=f,fill="white",anchor="lm")
        else:
            d.text((x+20,y+20),it["title"],font=f,fill=COLORS["ink"])
            body=ImageFont.truetype(FONT_PATH,it["fs"])
            lineheight=it["fs"]+11
            for j,line in enumerate(it["lines"]):
                d.text((x+22,y+62+j*lineheight),line,font=body,fill=COLORS["ink"])
    img.save(STEM.with_suffix(".png"), optimize=True)

if __name__=="__main__":
    drawio(); preview()
    print(STEM.with_suffix(".drawio"))
    print(STEM.with_suffix(".png"))
