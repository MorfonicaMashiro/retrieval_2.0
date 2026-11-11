"use strict";

const fs = require("fs");
const path = require("path");

const moduleRoot = process.argv[2];
const outputDir = process.argv[3];
if (!moduleRoot || !outputDir) {
  throw new Error("Usage: node render_forward_equations.js <node_modules> <output_dir>");
}

const sharp = require(path.join(moduleRoot, "sharp"));
const { mathjax } = require(path.join(moduleRoot, "mathjax-full/js/mathjax.js"));
const { TeX } = require(path.join(moduleRoot, "mathjax-full/js/input/tex.js"));
const { SVG } = require(path.join(moduleRoot, "mathjax-full/js/output/svg.js"));
const { liteAdaptor } = require(path.join(moduleRoot, "mathjax-full/js/adaptors/liteAdaptor.js"));
const { RegisterHTMLHandler } = require(path.join(moduleRoot, "mathjax-full/js/handlers/html.js"));
const { AllPackages } = require(path.join(moduleRoot, "mathjax-full/js/input/tex/AllPackages.js"));

const adaptor = liteAdaptor();
RegisterHTMLHandler(adaptor);
const tex = new TeX({ packages: AllPackages });
const svgOutput = new SVG({ fontCache: "local" });
const document = mathjax.document("", { InputJax: tex, OutputJax: svgOutput });

const equations = {
  model_input: String.raw`T_B=T_B(h_i,h_s,T_{\mathrm{sfc}},\mathrm{ice\ type},T_{\mathrm{water}},S_{\mathrm{water}})`,
  conductivity: String.raw`k_{\mathrm{ice}}=2.034+0.13\,\frac{S_{\mathrm{ice}}}{T_{\mathrm{ice}}-273.15},\qquad k_{\mathrm{snow}}=0.31\ \mathrm{W\,m^{-1}\,K^{-1}}`,
  heat_profile: String.raw`R_{\mathrm{th}}=\frac{h_s}{k_s}+\frac{h_i}{k_i},\quad q=\frac{T_w-T_{\mathrm{sfc}}}{R_{\mathrm{th}}},\quad T(z)=T_0+\frac{qz}{k}`,
  fyi_salinity: String.raw`S_{\mathrm{FYI}}=6.08e^{-5.81h_i}+7.409e^{-0.5228h_i},\qquad S_{\mathrm{FYI}}\ge 1.5\ \mathrm{ppt}`,
  myi_salinity: String.raw`z_k=\frac{k-1/2}{N},\qquad S_{i,k}=\frac{S_{\max}}{2}\left[1-\cos\!\left(\pi z_k^{\,a/(z_k+b)}\right)\right]`,
  myi_constants: String.raw`S_{\max}=3.2\ \mathrm{ppt},\qquad a=0.407,\qquad b=0.573`,
  snow_eps: String.raw`\begin{aligned}\varepsilon_s'&=1+0.7\rho_s+0.7\rho_s^2,\\\varepsilon_s''&=1.59\times10^6(0.52\rho_s+0.62\rho_s^2)\left(f^{-1}+1.23\times10^{-14}\sqrt f\right)e^{0.036T_s},\\\varepsilon_s&=\varepsilon_s'+i\varepsilon_s''\end{aligned}`,
  ice_eps: String.raw`\varepsilon_i=a_1+a_2V_b+i(a_3+a_4V_b)`,
  brine: String.raw`\rho_i=0.917-1.403\times10^{-4}T_i,\qquad V_b=1000\,\frac{\rho_iS_i}{F_1(T_i)-\rho_iS_iF_2(T_i)}`,
  water_eps: String.raw`\varepsilon_w=\varepsilon_\infty+\frac{\varepsilon_s(T,S)-\varepsilon_\infty}{1+i\omega\tau(T,S)}-i\frac{\sigma(T,S)}{\omega\varepsilon_0},\qquad \varepsilon_\infty=4.9`,
  fresnel: String.raw`q_j=\sqrt{\varepsilon_j-\sin^2\theta_0},\quad r_H=\frac{q_1-q_2}{q_1+q_2},\quad r_V=\frac{\varepsilon_2q_1-\varepsilon_1q_2}{\varepsilon_2q_1+\varepsilon_1q_2},\quad R=|r|^2`,
  attenuation: String.raw`\tau_j=\exp\!\left[-2k_0\left|\operatorname{Im}\sqrt{\varepsilon_j-\sin^2\theta_0}\right|d_j\right],\qquad k_0=\frac{2\pi f}{c_0}`,
  recursion: String.raw`\begin{aligned}A_j&=(1-\tau_j)(1+R_{j,b}\tau_j),\\U_j&=T_jA_j+(1-R_{j,b})U_{j+1}\tau_j,\\E_j&=A_j+(1-R_{j,b})E_{j+1}\tau_j\end{aligned}`,
  top_surface: String.raw`T_{B,\mathrm{surf}}=(1-R_{\mathrm{top}})U_1,\qquad e_{\mathrm{eff}}=(1-R_{\mathrm{top}})E_1`,
  toa: String.raw`T_{B,\mathrm{TOA}}=T_{B,\mathrm{surf}}+(1-e_{\mathrm{eff}})T_{B,\mathrm{cosm}}+\Delta T_{B,\mathrm{atm}}+\Delta T_{B,\mathrm{offset}}`,
  angular_mean: String.raw`\overline{T_B}=\frac{1}{2N_\theta}\sum_{p\in\{H,V\}}\sum_{\theta=0^\circ}^{40^\circ}T_B(p,\theta)`
};

function latexToSvg(latex) {
  const html = adaptor.outerHTML(document.convert(latex, { display: true }));
  const start = html.indexOf("<svg");
  const end = html.indexOf("</svg>");
  let svg = start >= 0 && end >= 0 ? html.slice(start, end + 6) : html;
  svg = svg.replace(/<\?xml[^>]*>/g, "");
  if (!/xmlns="http:\/\/www\.w3\.org\/2000\/svg"/.test(svg)) {
    svg = svg.replace(/<svg /, '<svg xmlns="http://www.w3.org/2000/svg" ');
  }
  svg = svg.replace(/currentColor/g, "#000000");
  return svg;
}

async function main() {
  fs.mkdirSync(outputDir, { recursive: true });
  for (const [name, latex] of Object.entries(equations)) {
    const svg = Buffer.from(latexToSvg(latex));
    await sharp(svg, { density: 300 })
      .trim({ background: "#ffffff" })
      .png()
      .toFile(path.join(outputDir, `${name}.png`));
  }
}

main().catch((error) => {
  process.stderr.write(String(error.stack || error));
  process.exit(1);
});
