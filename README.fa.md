# کتابخانه پژوهشی aoptics

**[English](README.md) | [فارسی](README.fa.md) | [Türkçe](README.tr.md)**

[![PyPI version](https://img.shields.io/pypi/v/aoptics.svg)](https://pypi.org/project/aoptics/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/aoptics.svg)](https://pypi.org/project/aoptics/)

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
- `aoptics.beams` — پرتوهای گاوسی، لاگر-گاوسی، هرمیت-گاوسی، بسل
- `aoptics.propagation` — طیف زاویه‌ای، فرنل، فرانهوفر، روش انتشار پرتو (BPM)، انتشار پرتو گاوسی با ماتریس ABCD
- `aoptics.interference` — تداخل دو/چندپرتویی، فابری-پرو، مایکلسون/مخ-زندر، دوشکاف، تداخل‌سنجی تغییر فاز، اسپکل، هولوگرافی
- `aoptics.diffraction` — توری‌ها، فرش تالبوت، الگوهای موآره، صفحات زون فرنل
- `aoptics.wavefront` — چندجمله‌ای‌های زرنیک، نسبت استرل، شبیه‌سازی شک-هارتمن
- `aoptics.quantum` — حالت‌های فاک/همدوس/فشرده/گرمایی/گربه، تابع ویگنر، آمار فوتون، حالت‌های بل

### فاز ۱ — اپتیک هندسی و الکترومغناطیسی (بورن-ولف فصل ۱،۳-۵،۱۴؛ صالح-تیچ فصل ۱،۶)
- `aoptics.geometrical.ray_tracing` — اپتیک ماتریسی پارکسیال، نقاط اصلی، ردیابی پرتوی واقعی سه‌بعدی دقیق
- `aoptics.geometrical.aberrations` — نظریه آبراسیون سایدل مرتبه‌سوم، ضرایب موج‌جبهه
- `aoptics.electromagnetic.fresnel` — ضرایب فرنل، زاویه بروستر/بحرانی
- `aoptics.electromagnetic.polarization` — حساب جونز، استوکس-مولر، کره پوانکاره
- `aoptics.electromagnetic.crystal_optics` — بیضی‌وار شاخص، دوشکستی، طراحی تیغه‌موج

### فاز ۲ — همدوسی و اپتیک آماری (بورن-ولف فصل ۱۰؛ گودمن)
- `aoptics.coherence.temporal` — قضیه وینر-خینچین، زمان/طول همدوسی
- `aoptics.coherence.spatial` — قضیه ون‌سیتریت-زرنیک، شعاع همدوسی
- `aoptics.coherence.speckle` — آمار اسپکل، شبیه‌سازی مونت‌کارلو

### فاز ۳ — فوتونیک (صالح-تیچ فصل ۷-۹،۲۱)
- `aoptics.photonics.waveguides` — مودهای TE موجبر لایه‌ای
- `aoptics.photonics.fibers` — عدد V، روزنه عددی، پاشندگی سلمایر
- `aoptics.photonics.resonators` — فابری-پرو، مود گاوسی خودسازگار
- `aoptics.photonics.nonlinear_optics` — تطبیق فاز SHG، اثر کر
- `aoptics.photonics.photonic_crystals` — روش ماتریس انتقال، آینه براگ

### فاز ۴ — کوانتوم اپتیک پیشرفته (گری-نایت، کامل)
- `aoptics.qoptics.operators` — عملگرهای میدان، جابجایی/فشرده‌سازی
- `aoptics.qoptics.jaynes_cummings` — مدل کامل جینز-کامینگز، فروپاشی-احیا
- `aoptics.qoptics.master_equation` — معادله مستر لیندبلاد
- `aoptics.qoptics.beamsplitter` — اثر هونگ-او-مندل
- `aoptics.qoptics.coherence_functions` — g²(۰)
- `aoptics.qoptics.entanglement` — آنتروپی فون‌نویمان، وفاداری
- `aoptics.qoptics.quantum_information` — تله‌پورت کوانتومی، BB84

### فاز ۵ — فوریه اپتیک، منابع/آشکارسازها، مدولاسیون، پراکندگی
- `aoptics.fourier_optics` — تابع مردمک، PSF، OTF/MTF، معیارهای تفکیک
- `aoptics.sources_detectors` — معادلات نرخ لیزر، نویز آشکارساز
- `aoptics.modulation` — الکترواپتیک (پاکلز)، آکوستواپتیک (براگ)
- `aoptics.scattering` — پراکندگی رایلی، اپتیک فلزات

### فاز ۶ — اپتیک غیرخطی طبق فصل‌های کتاب بوید
- `aoptics.nonlinear_optics_boyd.susceptibility` — بسط قطبش، قاعده میلر
- `aoptics.nonlinear_optics_boyd.coupled_wave_mixing` — اختلاط سه‌موج، روابط منلی-رو
- `aoptics.nonlinear_optics_boyd.two_level_atom` — اشباع اتم دوترازه
- `aoptics.nonlinear_optics_boyd.self_action` — خودکانونی‌سازی، فرمول ماربرگر
- `aoptics.nonlinear_optics_boyd.stimulated_scattering` — پراکندگی رامان/بریلوئن برانگیخته
- `aoptics.nonlinear_optics_boyd.multiphoton_absorption` — جذب دوفوتونی

## محدودیت‌های عمدی و آگاهانه دامنه

برخی موضوعات به‌عمد کنار گذاشته شدند چون فرمول‌های موجود بین منابع مختلف ریسک واقعی
تناقض قراردادی دارند (بهتر است حذف شود تا این‌که با حدس اشتباه ارائه شود): مدل گاوسی-شل،
نظریه کامل پراکندگی می (فقط حد رایلی پوشش داده شده)، دینامیک کامل قفل‌مد پالس فوق‌سریع،
و مؤلفه‌های کامل (u,v) معادلات بلاخ نوری.

## نصب
```bash
pip install aoptics
```

## راهنمای نگهداری: انتشار در PyPI

به فایل `docs/PUBLISHING.fa.md` مراجعه کنید — راهنمای کامل ساخت بسته، آپلود در PyPI/TestPyPI،
و راه‌اندازی انتشار خودکار از طریق GitHub Actions به‌طوری‌که هر انتشار (Release) جدید در گیت‌هاب
به‌صورت خودکار در PyPI منتشر شود.

## شروع سریع
```python
import aoptics as op
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
from aoptics.geometrical import ray_tracing as rt, aberrations as ab
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


