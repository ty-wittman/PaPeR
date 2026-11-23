// Print a rendered HTML file to PDF with Chromium: node print.js in.html out.pdf
const { chromium } = require('playwright');
(async () => {
  const [inp, out] = process.argv.slice(2);
  const browser = await chromium.launch();
  const page = await browser.newPage();
  await page.goto('file://' + require('path').resolve(inp), { waitUntil: 'networkidle' });
  await page.waitForFunction(() => document.querySelectorAll('.katex').length >= 0 && document.fonts.status === 'loaded');
  await page.waitForTimeout(400);
  await page.pdf({ path: out, format: 'Letter', printBackground: true, preferCSSPageSize: true,
    displayHeaderFooter: true, headerTemplate: '<span></span>',
    footerTemplate: '<div style="width:100%;text-align:center;font-family:serif;font-size:8.5pt;color:#555;"><span class="pageNumber"></span></div>',
    margin: { top: '22mm', bottom: '24mm', left: '22mm', right: '22mm' } });
  await browser.close();
})();
