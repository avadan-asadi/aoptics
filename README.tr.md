# aoptics Araştırma Kütüphanesi

**[English](README.md) | [فارسی](README.fa.md) | [Türkçe](README.tr.md)**

[![PyPI version](https://img.shields.io/pypi/v/aoptics.svg)](https://pypi.org/project/aoptics/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python](https://img.shields.io/pypi/pyversions/aoptics.svg)](https://pypi.org/project/aoptics/)

Klasik ve kuantum optik araştırma ve modellemesi için kapsamlı bir Python kütüphanesi.
Aşağıdaki temel kaynak kitapların bölümlerini sistematik biçimde kapsayacak şekilde
geliştirilmiştir:
- Born & Wolf, *Principles of Optics*
- Saleh & Teich, *Fundamentals of Photonics*
- Goodman, *Statistical/Fourier Optics*
- Gerry & Knight, *Introductory Quantum Optics*
- Boyd, *Nonlinear Optics*

**190 birim testi** kütüphane genelindeki fiziği doğrular (ders kitabı kapalı-form
sonuçlarıyla analitik karşılaştırmalar, enerji/olasılık korunumu, Monte Carlo
simülasyonlarının teoriyle uyumu ve aynı büyüklüğün bağımsız türetimleri arasında
çapraz doğrulama). Her modülün karşılık gelen bir `tests/test_*.py` dosyası vardır ve
her örnek betik uçtan uca çalıştırılabilir, bir grafik üretir.

## Modül haritası

### Temel çekirdek (dalga optiği / skaler kırınım / kuantum durumları)
- `aoptics.beams` — Gauss, Laguerre-Gauss, Hermite-Gauss, Bessel demetleri
- `aoptics.propagation` — açısal spektrum, Fresnel, Fraunhofer, demet yayılım yöntemi (BPM)
- `aoptics.interference` — çift/çok demetli girişim, Fabry-Perot, Michelson/Mach-Zehnder, holografi
- `aoptics.diffraction` — kırınım ağları, Talbot halıları, Moire desenleri, Fresnel bölge plakaları
- `aoptics.wavefront` — Zernike polinomları, Strehl oranı, Shack-Hartmann simülasyonu
- `aoptics.quantum` — Fock/uyumlu/sıkıştırılmış/termal/kedi durumları, Wigner fonksiyonu

### Faz 1 — Geometrik ve elektromanyetik optik (Born & Wolf Böl. 1,3-5,14; Saleh & Teich Böl. 1,6)
- `aoptics.geometrical.ray_tracing` — paraksiyel ABCD matris optiği, tam 3B ışın izleme
- `aoptics.geometrical.aberrations` — üçüncü dereceden Seidel aberasyon teorisi
- `aoptics.electromagnetic.fresnel` — Fresnel katsayıları, Brewster/kritik açı
- `aoptics.electromagnetic.polarization` — Jones hesabı, Stokes/Mueller hesabı, Poincare küresi
- `aoptics.electromagnetic.crystal_optics` — indeks elipsoidi, çift kırılım

### Faz 2 — Uyumluluk ve istatistiksel optik (Born & Wolf Böl. 10; Goodman)
- `aoptics.coherence.temporal` — Wiener-Khinchin teoremi, uyumluluk süresi/uzunluğu
- `aoptics.coherence.spatial` — Van Cittert-Zernike teoremi
- `aoptics.coherence.speckle` — speckle istatistikleri, Monte Carlo simülasyonu

### Faz 3 — Fotonik (Saleh & Teich Böl. 7-9,21)
- `aoptics.photonics.waveguides` — düzlemsel dalga kılavuzu TE modları
- `aoptics.photonics.fibers` — V-sayısı, sayısal açıklık, Sellmeier dağılımı
- `aoptics.photonics.resonators` — Fabry-Perot, kendiliğinden tutarlı Gauss modu
- `aoptics.photonics.nonlinear_optics` — SHG faz eşleştirme, Kerr etkisi
- `aoptics.photonics.photonic_crystals` — çok katmanlı transfer matrisi yöntemi

### Faz 4 — İleri kuantum optik (Gerry & Knight, tam metin)
- `aoptics.qoptics.operators` — alan operatörleri, yer değiştirme/sıkıştırma
- `aoptics.qoptics.jaynes_cummings` — tam Jaynes-Cummings modeli, çöküş-yeniden doğuş
- `aoptics.qoptics.master_equation` — Lindblad ana denklemi
- `aoptics.qoptics.beamsplitter` — Hong-Ou-Mandel etkisi
- `aoptics.qoptics.entanglement` — von Neumann entropisi, sadakat
- `aoptics.qoptics.quantum_information` — kuantum ışınlanma, BB84

### Faz 5 — Fourier optiği, kaynaklar/dedektörler, modülasyon, saçılma
- `aoptics.fourier_optics` — açıklık fonksiyonu, PSF, OTF/MTF, çözünürlük kriterleri
- `aoptics.sources_detectors` — lazer hız denklemleri, fotodedektör gürültüsü
- `aoptics.modulation` — elektro-optik (Pockels), akusto-optik (Bragg)
- `aoptics.scattering` — Rayleigh saçılması, metal optiği

### Faz 6 — Boyd'un kitap bölümlerine göre doğrusal olmayan optik
- `aoptics.nonlinear_optics_boyd.susceptibility` — polarizasyon açılımı, Miller kuralı
- `aoptics.nonlinear_optics_boyd.coupled_wave_mixing` — üç dalga karışımı, Manley-Rowe ilişkileri
- `aoptics.nonlinear_optics_boyd.two_level_atom` — iki seviyeli atom doygunluğu
- `aoptics.nonlinear_optics_boyd.self_action` — öz-odaklanma, Marburger formülü
- `aoptics.nonlinear_optics_boyd.stimulated_scattering` — uyarılmış Raman/Brillouin saçılması
- `aoptics.nonlinear_optics_boyd.multiphoton_absorption` — iki foton soğurma

## Bilinçli olarak dışarıda bırakılan konular

Bazı konular, kaynaklar arasında gerçek bir gösterim-kuralı riski taşıdığı için kasıtlı
olarak dışarıda bırakıldı (yanlış bir tahminle sunmaktansa atlamak daha iyidir): Gauss-Schell
modeli kısmi uyumlu demet yayılımı, tam Mie saçılma teorisi (yalnızca küçük parçacık
Rayleigh limiti kapsanmıştır), ayrıntılı mod kilitli ultra hızlı darbe dinamikleri.

## Kurulum
```bash
pip install aoptics
```

## Bakımcılar için: PyPI'ye yayınlama

Tam kılavuz için `docs/PUBLISHING.tr.md` dosyasına bakın: paketi oluşturma, PyPI/TestPyPI'ye
yükleme ve GitHub Actions üzerinden otomatik yayınlama kurulumu — böylece GitHub'da
oluşturulan her yeni Release otomatik olarak PyPI'ye yayınlanır.

## Hızlı başlangıç
```python
import aoptics as op
import numpy as np

# Laguerre-Gauss demeti oluştur
field = op.beams.laguerre_gaussian(p=0, l=2)
op.utils.show_field(field, title='LG(0,2) beam')

# Girişim deseni
lg = op.beams.laguerre_gaussian(0, 1)
ref = op.beams.gaussian_beam()
I, phase, vis = op.interference.two_beam_interference(lg, ref)
op.utils.show_interference(I, phase)

# Seidel aberasyonlu basit bir mercek (Faz 1)
from aoptics.geometrical import ray_tracing as rt, aberrations as ab
surfaces = [rt.Surface(radius=51.5, thickness=5.3, index=1.5168, semi_diameter=12.5),
            rt.Surface(radius=np.inf, thickness=0.0, index=1.0, semi_diameter=12.5)]
system = rt.OpticalSystem(surfaces)
print(system.cardinal_points())
```

## Geliştirme ve test
```bash
pip install -e .[dev]
pytest tests/ -q
python examples/phase6_boyd_nonlinear_demo.py
```

## Lisans

MIT — bkz. [LICENSE](LICENSE).


