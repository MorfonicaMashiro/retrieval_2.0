"""生成三张可编辑 draw.io 流程图与 2 倍分辨率 PPT 预览图。"""
from pathlib import Path
import xml.etree.ElementTree as ET
import json
from PIL import Image, ImageDraw, ImageFont

OUT = Path(__file__).resolve().parent
W, H, SCALE = 1200, 510, 2
FONT = r"C:\Windows\Fonts\msyh.ttc"
WHITE = '#FFFFFF'
INK = '#163428'
GREEN = '#236431'
GREEN_LIGHT = '#EAF5E9'
BLUE = '#EDF5FA'
BLUE_LINE = '#5286A3'
ORANGE = '#FFF3E5'
ORANGE_LINE = '#C78948'
PURPLE = '#F3EFFA'
PURPLE_LINE = '#8C75AB'
GRAY = '#F4F6F5'
GRAY_LINE = '#8B9D93'

def vertex(name, x, y, w, h, title, lines, fill=WHITE, line=GREEN,
           fs=18, hfs=20, left=True):
    return dict(name=name, x=x, y=y, w=w, h=h, title=title, lines=lines,
                fill=fill, line=line, fs=fs, hfs=hfs, left=left)

def connector(name, points, color=GREEN, dashed=False, arrow=True):
    return dict(name=name, points=points, color=color, dashed=dashed, arrow=arrow)

figures = [
    dict(stem='01_数据获取与配准', boxes=[
        vertex('smos',30,20,530,95,'SMOS 逐日亮温',[
            'ICDC L3B v3，源网格 896×608，EPSG:3413',
            '读取同日 TB、nPair；−999 → NaN'],BLUE,BLUE_LINE,18),
        vertex('oib',640,20,530,95,'OIB 逐航次数据',[
            '2012/13 IDCSI4；2014/15 Quick Look',
            '−99999 → NaN；FBice = FBsnow − hs'],BLUE,BLUE_LINE,18),
        vertex('smos_ease',30,150,530,95,'SMOS 最近邻重采样到原版 EASE',[
            '目标 1441×1441，格距 12,533.76 m',
            '目标中心逆映射至源像元；nPair≤0 保留缺测'],GREEN_LIGHT,GREEN,18),
        vertex('oib_ease',640,150,530,95,'OIB 航迹点投影与落格',[
            '经纬度 → EASE (x,y)，只变位置、不插值',
            'col=floor((x−xMin)/Δ)+1；row=floor((yMax−y)/Δ)+1'],GREEN_LIGHT,GREEN,18),
        vertex('match',135,280,930,105,'同日 3×3 窗口配准，得到一条中心格输入',[
            '以有效 TB 中心格汇集本格及周围 8 格 OIB 原始点；只保留单一 FYI/MYI、无开水/薄冰/灰冰的样本',
            '同一批 M 点分别取均值：FBsnow、hs、FBice、hi；输入 TBobs、冰型、两种干舷；hi/hs 留作验证'],ORANGE,ORANGE_LINE,17),
        vertex('counts',30,410,1140,95,'2012—2015 配对轨迹与网格数',[
            '2012：14 航次 / 13 有效 / 3,690 格     2013：10 / 9 / 4,389 格     2014：15 / 15 / 5,433 格     2015：9 / 9 / 2,852 格',
            '合计 48 航次、46 条有有效配对、16,364 格；论文年份 2012+2015：23 航次、22 条有效、6,542 格'],PURPLE,PURPLE_LINE,17),
    ], edges=[
        connector('s1',[(295,116),(295,149)]),connector('o1',[(905,116),(905,149)]),
        connector('s2',[(295,246),(295,262),(600,262),(600,279)]),
        connector('o2',[(905,246),(905,262),(600,262)],arrow=False),
        connector('match_count',[(600,386),(600,409)]),
    ]),
    dict(stem='02_多层前向亮温',boxes=[
        vertex('profile',30,20,545,184,'① 温度与盐度廓线',[
            '20 冰层 + 1 雪层；各层取中点温度',
            'Rth=hs/ks+hi/ki；q=(Tw−Ts)/Rth',
            'FYI：S=max[1.5,6.08e^(−5.81hi)+7.409e^(−0.5228hi)]',
            'MYI：S(z)=1.6[1−cos(π z^(0.407/(z+0.573)))]'],GREEN_LIGHT,GREEN,17),
        vertex('dielectric',625,20,545,184,'② 卤水体积与复介电常数',[
            'Vb=1000·ρpure(T)·S/F1(T)（‰）；基线不截断',
            'Vant FYI：εice=3.10+0.0084Vb+i(0.037+0.00445Vb)',
            'Vant MYI：εice=3.10+0.0084Vb+i(0.003+0.00435Vb)',
            '干雪 Tiuri；海水 Klein–Swift Debye'],GREEN_LIGHT,GREEN,17),
        vertex('optics',30,239,545,174,'③ 层界面与光学衰减',[
            'qj=√(εj−sin²θ)；Fresnel 反射 Rj=|rH/V|²',
            'τj=exp[−2k0|Im(qj)|dj]，k0=2πf/c',
            '同材数值层界面 R=0；物理界面用中点 ε',
            '单位温度源重复递推，得到有效发射率 e_eff'],BLUE,BLUE_LINE,18),
        vertex('rt',625,239,545,174,'④ Burke 自下而上亮温递推',[
            '海水边界 Ubottom=Tw(K)，向上穿过冰和雪',
            'Aj=(1−τj)(1+Rjτj)',
            'Uj=Tj(K)Aj+(1−Rj)τjU(j+1)',
            'TBθ,p=(1−Rtop)Utop+(1−e_eff)×2.7 K'],BLUE,BLUE_LINE,18),
        vertex('mean',30,445,1140,65,'输出及物理敏感性分支',[
            '0–40°、H/V 算术平均；实验：Vb 70‰ 截断、冰层 1/20/100、界面中点/边界、Burke/上层反射、雪分层'],PURPLE,PURPLE_LINE,17,19),
    ],edges=[
        connector('p_d',[(576,112),(624,112)]),
        connector('p_o',[(300,205),(300,238)]),
        connector('d_r',[(898,205),(898,238)]),
        connector('o_r',[(576,326),(624,326)]),
        connector('r_mean',[(898,414),(898,430),(600,430),(600,444)]),
    ]),
    dict(stem='03_反演评定与实验分支',boxes=[
        vertex('ice',30,20,545,152,'路线 A：TB + 冰面干舷 FBice',[
            '候选 hs：0–1 m，每 4 cm 扫描',
            'hi=(ρw·FBice+ρs·hs)/(ρw−ρi)',
            'g(hs)=F(hi,hs,type)−TBobs；保留物理解',
            '扫描命中或变号区间进入二分'],ORANGE,ORANGE_LINE,18),
        vertex('snow',625,20,545,152,'路线 B：TB + 雪面干舷 FBsnow',[
            '候选 hs：0–min(1 m,FBsnow)，每 1 cm 全扫',
            'hi=[ρw·FBsnow−(ρw−ρs)hs]/(ρw−ρi)',
            'FBice=FBsnow−hs≥0；每个变号区间二分',
            '保留全部根，验证默认取较小 hs 根'],ORANGE,ORANGE_LINE,18),
        vertex('bisect',130,205,940,101,'共同求根与精度控制',[
            '|g|≤0.01 K 或雪深区间≤1e−5 m，最多 60 次；入表前再次校验 TB 残差',
            '保存无根/双根、hi、hs、TBsim；不提前舍入；每条同日配对独立求解'],BLUE,BLUE_LINE,18),
        vertex('metrics',30,340,545,135,'误差评定',[
            'OIB 平均 hi/hs 作参照；TBOIB=F(hiOIB,hsOIB)',
            '报告 n、bias、RMSE、r、R²、无根率/双根率',
            '全部、逐年、nOIB 前 5%、1.5 K 诊断子集',
            '1.5 K 子集不用于实际反演前筛选'],GREEN_LIGHT,GREEN,18),
        vertex('branches',625,340,545,129,'空间算子与后续测试分支',[
            'stride=1 重叠；stride=3 不重叠',
            'A: F(mean hi,mean hs)；B: mean F(逐点 hi,hs)',
            'Aw: 加权参数；Bw: 加权逐点 TB 后反演',
            '权重：等面积子网格或中心反距离 w=d^(−p)'],PURPLE,PURPLE_LINE,17),
    ],edges=[
        connector('ice_b',[(300,173),(300,188),(600,188),(600,204)]),
        connector('snow_b',[(898,173),(898,188),(600,188)],arrow=False),
        connector('b_m',[(600,307),(600,322),(302,322),(302,339)]),
        connector('b_exp',[(600,322),(898,322),(898,339)],color=PURPLE_LINE,dashed=True),
    ]),
]

def hexrgb(s): return tuple(int(s[i:i+2],16) for i in (1,3,5))

def render_png(spec, target):
    im = Image.new('RGB',(W*SCALE,H*SCALE),WHITE)
    d = ImageDraw.Draw(im)
    def line(coords, color, width=2):
        d.line([(int(x*SCALE),int(y*SCALE)) for x,y in coords],fill=color,width=width*SCALE)
    for e in spec['edges']:
        pts=e['points']
        for a,b in zip(pts,pts[1:]):
            if e['dashed']:
                n=max(1,int((abs(a[0]-b[0])+abs(a[1]-b[1]))/8))
                for i in range(0,n,2):
                    p=(a[0]+(b[0]-a[0])*i/n,a[1]+(b[1]-a[1])*i/n)
                    q=(a[0]+(b[0]-a[0])*min(i+1,n)/n,a[1]+(b[1]-a[1])*min(i+1,n)/n)
                    line([p,q],e['color'])
            else: line([a,b],e['color'])
        if e['arrow']:
            x,y=pts[-1]; ax,ay=pts[-2]
            if y>ay: tri=[(x,y),(x-6,y-10),(x+6,y-10)]
            elif x>ax: tri=[(x,y),(x-10,y-6),(x-10,y+6)]
            else: tri=[(x,y),(x+10,y-6),(x+10,y+6)]
            d.polygon([(int(u*SCALE),int(v*SCALE)) for u,v in tri],fill=e['color'])
    for b in spec['boxes']:
        x,y,w,h=(b[k]*SCALE for k in ('x','y','w','h'))
        d.rounded_rectangle((x,y,x+w,y+h),radius=9*SCALE,
                            fill=b['fill'],outline=b['line'],width=2*SCALE)
        heading=ImageFont.truetype(FONT,b['hfs']*SCALE)
        body=ImageFont.truetype(FONT,b['fs']*SCALE)
        d.text((x+15*SCALE,y+11*SCALE),b['title'],font=heading,fill=INK)
        for i,t in enumerate(b['lines']):
            yy=y+(42+i*(b['fs']+4))*SCALE
            d.text((x+16*SCALE,yy),t,font=body,fill=INK)
    im.save(target,optimize=True)

def render_drawio(spec,target):
    root=ET.Element('mxfile',{'host':'app.diagrams.net','agent':'Codex','version':'24.7.17','pages':'1'})
    diag=ET.SubElement(root,'diagram',{'id':spec['stem'],'name':spec['stem']})
    graph=ET.SubElement(diag,'mxGraphModel',{'dx':str(W),'dy':str(H),'grid':'0','page':'1','pageScale':'1','pageWidth':str(W),'pageHeight':str(H),'math':'0','shadow':'0'})
    cells=ET.SubElement(graph,'root')
    ET.SubElement(cells,'mxCell',{'id':'0'})
    ET.SubElement(cells,'mxCell',{'id':'1','parent':'0'})
    for b in spec['boxes']:
        val='<b>'+b['title']+'</b><br>'+'<br>'.join(b['lines'])
        sty=('rounded=1;arcSize=8;whiteSpace=wrap;html=1;'
             f'fillColor={b["fill"]};strokeColor={b["line"]};strokeWidth=2;'
             f'fontSize={b["fs"]};fontColor={INK};fontFamily=Microsoft YaHei;'
             'align=left;verticalAlign=middle;spacingLeft=14;spacingRight=8;')
        c=ET.SubElement(cells,'mxCell',{'id':b['name'],'value':val,'style':sty,'vertex':'1','parent':'1'})
        ET.SubElement(c,'mxGeometry',{'x':str(b['x']),'y':str(b['y']),'width':str(b['w']),'height':str(b['h']),'as':'geometry'})
    for e in spec['edges']:
        sty=('edgeStyle=orthogonalEdgeStyle;rounded=0;html=1;'
             f'strokeColor={e["color"]};strokeWidth=2;dashed={int(e["dashed"])};'
             f'endArrow={"block" if e["arrow"] else "none"};endFill=1;endSize=7;')
        c=ET.SubElement(cells,'mxCell',{'id':e['name'],'value':'','style':sty,'edge':'1','parent':'1'})
        g=ET.SubElement(c,'mxGeometry',{'relative':'1','as':'geometry'})
        ET.SubElement(g,'mxPoint',{'x':str(e['points'][0][0]),'y':str(e['points'][0][1]),'as':'sourcePoint'})
        ET.SubElement(g,'mxPoint',{'x':str(e['points'][-1][0]),'y':str(e['points'][-1][1]),'as':'targetPoint'})
        if len(e['points'])>2:
            arr=ET.SubElement(g,'Array',{'as':'points'})
            for x,y in e['points'][1:-1]: ET.SubElement(arr,'mxPoint',{'x':str(x),'y':str(y)})
    ET.indent(root,space='  ')
    ET.ElementTree(root).write(target,encoding='utf-8',xml_declaration=True)

def check_text(spec):
    problems=[]
    for b in spec['boxes']:
        hfont=ImageFont.truetype(FONT,b['hfs'])
        font=ImageFont.truetype(FONT,b['fs'])
        if hfont.getlength(b['title'])>b['w']-30: problems.append((b['name'],'title width'))
        for i,t in enumerate(b['lines']):
            if font.getlength(t)>b['w']-32: problems.append((b['name'],f'line {i+1} width'))
            if 42+i*(b['fs']+4)+b['fs']+5>b['h']: problems.append((b['name'],f'line {i+1} height'))
    return problems

if __name__=='__main__':
    OUT.mkdir(parents=True,exist_ok=True)
    with (OUT/'三页流程图_元素.json').open('w',encoding='utf-8') as f:
        json.dump(figures,f,ensure_ascii=False,indent=2)
    for spec in figures:
        errors=check_text(spec)
        if errors: raise RuntimeError(f"{spec['stem']} overflow: {errors}")
        render_png(spec,OUT/(spec['stem']+'.png'))
        render_drawio(spec,OUT/(spec['stem']+'.drawio'))
        print(spec['stem'],'OK')
