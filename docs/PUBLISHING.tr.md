# Yayınlama Kılavuzu (Türkçe)

Bu kılavuz, **opticspy-research** paketini oluşturmak, PyPI'ye yüklemek ve her yeni GitHub
Release'in kendiliğinden PyPI'ye yayınlanmasını sağlayan **otomatik yayınlama** sistemini
kurmak için gereken her şeyi kapsar -- ilk kurulumdan sonra manuel `twine upload` gerekmez.

`import opticspy` ile çağrılan Python paketi, PyPI'de **`opticspy-research`** adıyla
dağıtılır (sade `opticspy` adı zaten başka, ilgisiz bir proje tarafından alınmıştır).

---

## 0. Tek seferlik ön koşullar

1. Bir [PyPI](https://pypi.org/account/register/) hesabı (ve prova yüklemeleri için önerilen
   bir [TestPyPI](https://test.pypi.org/account/register/) hesabı).
2. Bu projenin içeriğini barındıran bir GitHub deposu (bu zip'in içeriğini yeni bir depoya
   push edin).
3. İlk gerçek yüklemeden önce **`pyproject.toml` dosyasını düzenleyin**: `[project.urls]`
   bölümündeki her `OWNER` ifadesini gerçek GitHub kullanıcı adınız/organizasyonunuzla
   değiştirin.
4. Python 3.9+ ve yerel olarak `pip install build twine` (zaten `pip install -e .[dev]` ile
   kurulur).

---

## 1. Paketi yerel olarak oluşturma ve kontrol etme

```bash
cd opticspy-research/          # pyproject.toml dosyasının bulunduğu klasör
python -m pip install --upgrade build twine
python -m build                # dist/*.whl ve dist/*.tar.gz dosyalarını oluşturur
twine check dist/*             # meta veri/uzun açıklama render kontrolü yapar
```

Önce test paketini çalıştırarak her şeyin doğru çalıştığından emin olun:
```bash
pip install -e .[dev]
pytest tests/ -q
```

---

## 2. İlk yükleme: manuel olarak, önce TestPyPI sonra PyPI

Gerçek yayından önce her şeyin doğru render edildiğinden ve kurulduğundan emin olmak için
önce TestPyPI'ye yüklemek iyi bir uygulamadır.

```bash
# 1) TestPyPI'ye yükle
twine upload --repository testpypi dist/*
# (istendiğinde TestPyPI API token'ınızı girin, kullanıcı adı: __token__)

# 2) Doğrulamak için temiz bir venv'de TestPyPI'den kurun
python -m venv /tmp/test_env && source /tmp/test_env/bin/activate
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ opticspy-research
python -c "import opticspy; print(opticspy.__version__)"
deactivate

# 3) Her şey doğruysa, gerçek PyPI'ye yükleyin
twine upload dist/*
# (kullanıcı adı: __token__, şifre: PyPI API token'ınız)
```

API token alma: PyPI hesap ayarları -> "API tokens" -> "Add API token". İlk yükleme hesap
genelinde bir token gerektirir; sonrasında token kapsamını yalnızca bu projeyle
sınırlayabilirsiniz.

Bu ilk manuel yüklemeden sonra `pip install opticspy-research` herkes için çalışır.

---

## 3. Her gelecek sürüm için otomatik yayınlama (GitHub Actions)

Bu depo zaten `.github/workflows/publish.yml` dosyasını içerir. Bir **GitHub Release**
yayınladığınızda otomatik olarak çalışır, paketi oluşturur ve PyPI'nin önerdiği güncel
yöntem olan **Trusted Publishing** (OpenID Connect) ile -- GitHub'da hiçbir gizli API
token'ı saklamaya gerek kalmadan -- PyPI'ye yükler.

### 3.1 Tek seferlik kurulum: GitHub'ı PyPI'ye bağlama (Trusted Publishing)

1. PyPI'de projenize gidin:
   `https://pypi.org/manage/project/opticspy-research/settings/publishing/`
   (projenin var olması için Bölüm 2'deki ilk manuel yüklemeyi yapmış olmanız gerekir).
2. "Trusted Publishers" altında "Add a new publisher"a tıklayın ve şunları doldurun:
   - Owner: GitHub kullanıcı adınız/organizasyonunuz
   - Repository name: depo adınız (örn. `opticspy-research`)
   - Workflow name: `publish.yml`
   - Environment name: `pypi`
3. Kaydedin. Bu kadar -- hiçbir yere token kopyalamaya gerek yok.

### 3.2 Bir release yayınının otomatik tetiklemesi

1. **`opticspy/__init__.py`** dosyasındaki sürüm numarasını artırın (`__version__ = "2.2.0"`
   gibi) -- bu tek gerçek kaynaktır; `pyproject.toml` bunu otomatik olarak okur.
2. **`CHANGELOG.md`** dosyasının başına ilgili bir girdi ekleyin.
3. Commit edin, sonra etiketleyip push edin:
   ```bash
   git add -A && git commit -m "Release v2.2.0"
   git tag v2.2.0
   git push && git push --tags
   ```
4. GitHub'da "Releases" -> "Draft a new release"e gidin, az önce push ettiğiniz etiketi
   (örn. `v2.2.0`) seçin, sürüm notlarını yazın (veya CHANGELOG.md'den kopyalayın) ve
   "Publish release"e tıklayın.
5. Release'i yayınlamak otomatik olarak `.github/workflows/publish.yml` dosyasını tetikler;
   bu da paketi oluşturur ve bir-iki dakika içinde PyPI'ye yükler. Başka bir işlem gerekmez.

### 3.3 Alternatif: API token yöntemi (Trusted Publishing kullanmak istemiyorsanız)

OIDC yerine klasik bir token kullanmayı tercih ediyorsanız:
1. Bu projeyle sınırlı bir PyPI API token oluşturun.
2. GitHub deponuzda: Settings -> Secrets and variables -> Actions -> "New repository
   secret", adını `PYPI_API_TOKEN` yapın, token'ı yapıştırın.
3. `.github/workflows/publish.yml` içindeki `publish` job adımlarını şununla değiştirin:
   ```yaml
   - uses: pypa/gh-action-pypi-publish@release/v1
     with:
       password: ${{ secrets.PYPI_API_TOKEN }}
   ```
   (ve yalnızca OIDC için gereken `permissions: id-token: write` bloğunu kaldırın).

---

## 4. Sürekli entegrasyon (her push/PR'da testler)

`.github/workflows/tests.yml`, `main` dalına yapılan her push ve pull request'te testleri
Python 3.9-3.12 üzerinde çalıştırır, böylece regresyonlar bir sürüm etiketlenmeden önce
yakalanır.

---

## 5. Her yeni sürüm için hızlı kontrol listesi

- [ ] Kod değişiklikleri `main`e birleştirildi, tüm testler yeşil.
- [ ] `opticspy/__init__.py`: `__version__` artırıldı.
- [ ] `CHANGELOG.md`: yeni girdi eklendi.
- [ ] `git tag vX.Y.Z && git push --tags`
- [ ] GitHub'da o etiketten "Draft a new release" -> "Publish release".
- [ ] Actions sekmesinde `publish.yml` çalışmasının yeşile dönmesini izleyin.
- [ ] `pip install --upgrade opticspy-research` ile yeni sürümü doğrulayın.
