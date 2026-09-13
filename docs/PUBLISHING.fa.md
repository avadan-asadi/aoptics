# راهنمای انتشار (فارسی)

این راهنما تمام مراحل لازم برای ساخت بسته **opticspy-research**، آپلود آن در PyPI، و
راه‌اندازی **انتشار خودکار** را پوشش می‌دهد — به‌طوری‌که هر انتشار (Release) جدید در گیت‌هاب
خودش را به‌صورت خودکار در PyPI منتشر کند، بدون نیاز به `twine upload` دستی پس از راه‌اندازی اولیه.

پکیج پایتونی که با `import opticspy` فراخوانی می‌شود، در PyPI با نام **`opticspy-research`**
منتشر می‌شود (نام ساده `opticspy` قبلاً توسط پروژه دیگری اشغال شده است).

---

## ۰. پیش‌نیازهای یک‌باره

۱. یک حساب کاربری در [PyPI](https://pypi.org/account/register/) (و ترجیحاً یک حساب در
   [TestPyPI](https://test.pypi.org/account/register/) برای تست اولیه).
۲. یک مخزن (Repository) در گیت‌هاب که محتوای این zip در آن قرار گیرد.
۳. پیش از اولین آپلود واقعی، **فایل `pyproject.toml` را ویرایش کنید**: هر `OWNER` را در بخش
   `[project.urls]` با نام کاربری/سازمان واقعی گیت‌هاب خود جایگزین کنید.
۴. پایتون ۳.۹ به بالا و `pip install build twine` روی سیستم محلی (که با `pip install -e .[dev]` نصب می‌شود).

---

## ۱. ساخت و بررسی بسته به‌صورت محلی

```bash
cd opticspy-research/          # پوشه‌ای که pyproject.toml در آن است
python -m pip install --upgrade build twine
python -m build                # فایل‌های dist/*.whl و dist/*.tar.gz را می‌سازد
twine check dist/*             # صحت متادیتا/رندر توضیحات را بررسی می‌کند
```

قبل از هر چیز، مجموعه تست‌ها را اجرا کنید تا مطمئن شوید همه‌چیز درست کار می‌کند:
```bash
pip install -e .[dev]
pytest tests/ -q
```

---

## ۲. اولین آپلود: دستی، از طریق TestPyPI و سپس PyPI

بهتر است ابتدا در TestPyPI آپلود کنید تا از درستی رندر و نصب مطمئن شوید، سپس در PyPI واقعی منتشر کنید.

```bash
# ۱) آپلود در TestPyPI
twine upload --repository testpypi dist/*
# (نام کاربری: __token__، رمز عبور: توکن API خود در TestPyPI)

# ۲) نصب از TestPyPI در یک محیط مجازی تمیز برای تأیید
python -m venv /tmp/test_env && source /tmp/test_env/bin/activate
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ opticspy-research
python -c "import opticspy; print(opticspy.__version__)"
deactivate

# ۳) اگر همه‌چیز درست بود، در PyPI واقعی آپلود کنید
twine upload dist/*
# (نام کاربری: __token__، رمز عبور: توکن API خود در PyPI)
```

**دریافت توکن API:** تنظیمات حساب PyPI ← «API tokens» ← «Add API token». پس از اولین آپلود
(که نیاز به توکن در سطح کل حساب دارد)، می‌توانید دامنه توکن را محدود به همین پروژه کنید.

پس از این آپلود دستی اول، دستور `pip install opticspy-research` برای همه کار می‌کند.

---

## ۳. انتشار خودکار برای هر نسخه بعدی (GitHub Actions)

این مخزن از قبل شامل فایل `.github/workflows/publish.yml` است. این فایل هرگاه یک
**انتشار (Release) گیت‌هاب** منتشر کنید، به‌صورت خودکار اجرا می‌شود، بسته را می‌سازد و با
روش **Trusted Publishing** (پروتکل OpenID Connect) — روش توصیه‌شده فعلی PyPI که نیازی به
ذخیره‌ی هیچ توکن مخفی در گیت‌هاب ندارد — آن را در PyPI منتشر می‌کند.

### ۳.۱ راه‌اندازی یک‌باره: اتصال گیت‌هاب به PyPI (Trusted Publishing)

۱. به صفحه پروژه خود در PyPI بروید:
   `https://pypi.org/manage/project/opticspy-research/settings/publishing/`
   (باید آپلود دستی اول از بخش ۲ را قبلاً انجام داده باشید تا پروژه وجود داشته باشد).
۲. در بخش «Trusted Publishers»، روی «Add a new publisher» کلیک کنید و موارد زیر را وارد کنید:
   - Owner: نام کاربری/سازمان گیت‌هاب شما
   - Repository name: نام مخزن شما (مثلاً `opticspy-research`)
   - Workflow name: `publish.yml`
   - Environment name: `pypi`
۳. ذخیره کنید. همین — نیازی به کپی هیچ توکنی نیست.

### ۳.۲ نحوه فعال‌شدن انتشار با یک Release

۱. شماره نسخه را در فایل **`opticspy/__init__.py`** افزایش دهید (مثلاً `__version__ = "2.2.0"`) —
   این تنها منبع حقیقت نسخه است؛ `pyproject.toml` آن را به‌صورت خودکار می‌خواند.
۲. یک ورودی متناظر در ابتدای فایل **`CHANGELOG.md`** اضافه کنید.
۳. کامیت کنید، سپس تگ بزنید و پوش کنید:
   ```bash
   git add -A && git commit -m "Release v2.2.0"
   git tag v2.2.0
   git push && git push --tags
   ```
۴. در گیت‌هاب، به بخش «Releases» ← «Draft a new release» بروید، تگی که پوش کردید
   (مثلاً `v2.2.0`) را انتخاب کنید، توضیحات انتشار را بنویسید (یا از CHANGELOG.md کپی کنید)،
   و روی «Publish release» کلیک کنید.
۵. انتشار Release به‌صورت خودکار فایل `.github/workflows/publish.yml` را اجرا می‌کند که بسته
   را می‌سازد و ظرف یکی دو دقیقه در PyPI منتشر می‌کند. نیازی به کار دیگری نیست.

### ۳.۳ روش جایگزین: توکن API (اگر نمی‌خواهید از Trusted Publishing استفاده کنید)

اگر ترجیح می‌دهید از روش کلاسیک توکن به‌جای OIDC استفاده کنید:
۱. یک توکن API در PyPI محدود به این پروژه بسازید.
۲. در مخزن گیت‌هاب: Settings ← Secrets and variables ← Actions ← «New repository secret»،
   نام آن را `PYPI_API_TOKEN` بگذارید و توکن را وارد کنید.
۳. مراحل job با نام `publish` در `.github/workflows/publish.yml` را با این جایگزین کنید:
   ```yaml
   - uses: pypa/gh-action-pypi-publish@release/v1
     with:
       password: ${{ secrets.PYPI_API_TOKEN }}
   ```
   (و بخش `permissions: id-token: write` را که فقط برای OIDC لازم است حذف کنید).

---

## ۴. یکپارچه‌سازی پیوسته (تست در هر push/PR)

فایل `.github/workflows/tests.yml` در هر push و pull request به شاخه `main`، تست‌ها را روی
پایتون ۳.۹ تا ۳.۱۲ اجرا می‌کند تا مشکلات پیش از تگ‌زدن هر نسخه شناسایی شوند.

---

## ۵. چک‌لیست سریع برای هر انتشار جدید

- [ ] تغییرات کد در `main` ادغام شده و همه تست‌ها سبز است.
- [ ] `opticspy/__init__.py`: مقدار `__version__` افزایش یافته.
- [ ] `CHANGELOG.md`: ورودی جدید اضافه شده.
- [ ] `git tag vX.Y.Z && git push --tags`
- [ ] در گیت‌هاب «Draft a new release» از همان تگ ← «Publish release».
- [ ] در تب Actions منتظر بمانید تا اجرای `publish.yml` سبز شود.
- [ ] با `pip install --upgrade opticspy-research` نسخه جدید را تأیید کنید.
