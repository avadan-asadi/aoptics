# کتابخانه پژوهشی OpticsPy

**[English](README.md) | [فارسی](README.fa.md) | [Türkçe](README.tr.md)**

[![PyPI version](https://img.shields.io/pypi/v/opticspy-research.svg)](https://pypi.org/project/opticspy-research/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/opticspy-research.svg)](https://pypi.org/project/opticspy-research/)

یک کتابخانه جامع پایتون برای پژوهش و مدل‌سازی اپتیک کلاسیک و کوانتومی، که به‌صورت نظام‌مند
فصل‌های کتاب‌های زیر را پوشش می‌دهد:
- *Principles of Optics* نوشته بورن و ولف (Born & Wolf)
- *Fundamentals of Photonics* نوشته صالح و تیچ (Saleh & Teich)
- *Statistical/Fourier Optics* نوشته گودمن (Goodman)
- *Introductory Quantum Optics* نوشته گری و نایت (Gerry & Knight)
- *Nonlinear Optics* نوشته بوید (Boyd)

**۱۹۰ تست واحد** صحت فیزیکی کل کتابخانه را تأیید می‌کنند (مقایسه تحلیلی با نتایج بسته کتاب‌های
درسی، پایستگی انرژی/احتمال، تطابق شبیه‌سازی مونت‌کارلو با نظریه، و اعتبارسنجی متقابل بین
اشتقاق‌های مستقل یک کمیت). هر ماژول یک فایل `tests/test_*.py` متناظر دارد و هر اسکریپت مثال
به‌صورت کامل اجرا می‌شود و یک نمودار تولید می‌کند.

## نقشه ماژول‌ها

### هسته اصلی (اپتیک موج / پراش اسکالر / حالت‌های کوانتومی)
- `opticspy.beams` — پرتوهای گاوسی، لاگر-گاوسی، هرمیت-گاوسی، بسل
- `opticspy.propagation` — طیف زاویه‌ای، فرنل، فرانهوفر، روش انتشار پرتو (BPM)، انتشار پرتو گاوسی با ماتریس ABCD
- `opticspy.interference` — تداخل دو/چندپرتویی، فابری-پرو، مایکلسون/مخ-زندر، دوشکاف، تداخل‌سنجی تغییر فاز، اسپکل، هولوگرافی
- `opticspy.diffraction` — توری‌ها، فرش تالبوت، الگوهای موآره، صفحات زون فرنل
- `opticspy.wavefront` — چندجمله‌ای‌های زرنیک، نسبت استرل، شبیه‌سازی شک-هارتمن
- `opticspy.quantum` — حالت‌های فاک/همدوس/فشرده/گرمایی/گربه، تابع ویگنر، آمار فوتون، حالت‌های بل

### فاز ۱ — اپتیک هندسی و الکترومغناطیسی (بورن-ولف فصل ۱،۳-۵،۱۴؛ صالح-تیچ فصل ۱،۶)
- `opticspy.geometrical.ray_tracing` — اپتیک ماتریسی پارکسیال، نقاط اصلی، ردیابی پرتوی واقعی سه‌بعدی دقیق
- `opticspy.geometrical.aberrations` — نظریه آبراسیون سایدل مرتبه‌سوم، ضرایب موج‌جبهه
- `opticspy.electromagnetic.fresnel` — ضرایب فرنل، زاویه بروستر/بحرانی
- `opticspy.electromagnetic.polarization` — حساب جونز، استوکس-مولر، کره پوانکاره
- `opticspy.electromagnetic.crystal_optics` — بیضی‌وار شاخص، دوشکستی، طراحی تیغه‌موج

### فاز ۲ — همدوسی و اپتیک آماری (بورن-ولف فصل ۱۰؛ گودمن)
- `opticspy.coherence.temporal` — قضیه وینر-خینچین، زمان/طول همدوسی
- `opticspy.coherence.spatial` — قضیه ون‌سیتریت-زرنیک، شعاع همدوسی
- `opticspy.coherence.speckle` — آمار اسپکل، شبیه‌سازی مونت‌کارلو

### فاز ۳ — فوتونیک (صالح-تیچ فصل ۷-۹،۲۱)
- `opticspy.photonics.waveguides` — مودهای TE موجبر لایه‌ای
- `opticspy.photonics.fibers` — عدد V، روزنه عددی، پاشندگی سلمایر
- `opticspy.photonics.resonators` — فابری-پرو، مود گاوسی خودسازگار
- `opticspy.photonics.nonlinear_optics` — تطبیق فاز SHG، اثر کر
- `opticspy.photonics.photonic_crystals` — روش ماتریس انتقال، آینه براگ

### فاز ۴ — کوانتوم اپتیک پیشرفته (گری-نایت، کامل)
- `opticspy.qoptics.operators` — عملگرهای میدان، جابجایی/فشرده‌سازی
- `opticspy.qoptics.jaynes_cummings` — مدل کامل جینز-کامینگز، فروپاشی-احیا
- `opticspy.qoptics.master_equation` — معادله مستر لیندبلاد
- `opticspy.qoptics.beamsplitter` — اثر هونگ-او-مندل
- `opticspy.qoptics.coherence_functions` — g²(۰)
- `opticspy.qoptics.entanglement` — آنتروپی فون‌نویمان، وفاداری
- `opticspy.qoptics.quantum_information` — تله‌پورت کوانتومی، BB84

### فاز ۵ — فوریه اپتیک، منابع/آشکارسازها، مدولاسیون، پراکندگی
- `opticspy.fourier_optics` — تابع مردمک، PSF، OTF/MTF، معیارهای تفکیک
- `opticspy.sources_detectors` — معادلات نرخ لیزر، نویز آشکارساز
- `opticspy.modulation` — الکترواپتیک (پاکلز)، آکوستواپتیک (براگ)
- `opticspy.scattering` — پراکندگی رایلی، اپتیک فلزات

### فاز ۶ — اپتیک غیرخطی طبق فصل‌های کتاب بوید
- `opticspy.nonlinear_optics_boyd.susceptibility` — بسط قطبش، قاعده میلر
- `opticspy.nonlinear_optics_boyd.coupled_wave_mixing` — اختلاط سه‌موج، روابط منلی-رو
- `opticspy.nonlinear_optics_boyd.two_level_atom` — اشباع اتم دوترازه
- `opticspy.nonlinear_optics_boyd.self_action` — خودکانونی‌سازی، فرمول ماربرگر
- `opticspy.nonlinear_optics_boyd.stimulated_scattering` — پراکندگی رامان/بریلوئن برانگیخته
- `opticspy.nonlinear_optics_boyd.multiphoton_absorption` — جذب دوفوتونی

## محدودیت‌های عمدی و آگاهانه دامنه

برخی موضوعات به‌عمد کنار گذاشته شدند چون فرمول‌های موجود بین منابع مختلف ریسک واقعی
تناقض قراردادی دارند (بهتر است حذف شود تا این‌که با حدس اشتباه ارائه شود): مدل گاوسی-شل،
نظریه کامل پراکندگی می (فقط حد رایلی پوشش داده شده)، دینامیک کامل قفل‌مد پالس فوق‌سریع،
و مؤلفه‌های کامل (u,v) معادلات بلاخ نوری.

## نصب
```bash
pip install opticspy-research
```

## راهنمای نگهداری: انتشار در PyPI

به فایل `docs/PUBLISHING.fa.md` مراجعه کنید — راهنمای کامل ساخت بسته، آپلود در PyPI/TestPyPI،
و راه‌اندازی انتشار خودکار از طریق GitHub Actions به‌طوری‌که هر انتشار (Release) جدید در گیت‌هاب
به‌صورت خودکار در PyPI منتشر شود.

## شروع سریع
```python
import opticspy as op
import numpy as np

# تولید پرتو لاگر-گاوسی
field = op.beams.laguerre_gaussian(p=0, l=2)
op.utils.show_field(field, title='LG(0,2) beam')

# الگوی تداخل
lg = op.beams.laguerre_gaussian(0, 1)
ref = op.beams.gaussian_beam()
I, phase, vis = op.interference.two_beam_interference(lg, ref)
op.utils.show_interference(I, phase)

# یک لنز تک‌ساده با آبراسیون سایدل (فاز ۱)
from opticspy.geometrical import ray_tracing as rt, aberrations as ab
surfaces = [rt.Surface(radius=51.5, thickness=5.3, index=1.5168, semi_diameter=12.5),
            rt.Surface(radius=np.inf, thickness=0.0, index=1.0, semi_diameter=12.5)]
system = rt.OpticalSystem(surfaces)
print(system.cardinal_points())
```

## توسعه و تست
```bash
pip install -e .[dev]
pytest tests/ -q
python examples/phase6_boyd_nonlinear_demo.py
```

## لایسنس

MIT — به فایل [LICENSE](LICENSE) مراجعه کنید.
