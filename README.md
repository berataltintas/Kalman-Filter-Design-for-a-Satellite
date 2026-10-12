# Kalman Filter Design for a Satellite

Bu depo, bir uydunun tutum kestirimi için MATLAB tabanlı bir simülasyon ve çoklu **Multiplicative Extended Kalman Filter (MEKF)** yaklaşımı içerir.

## İçerik

Projede üç farklı kestirim akışı birlikte çalışır:

- **Coarse mode (Güneş + Manyetometre):** TRIAD ile başlangıç ve vektör ölçümlerinden MEKF güncellemesi.
- **Eclipse mode (Jiroskop + Manyetometre):** Güneş sensörü yokken (tutulma/gölge), jiroskop-bias ile manyetik ölçümden MEKF güncellemesi.
- **Fine mode (Yıldız İzleyici + Jiroskop):** Yüksek doğruluklu yıldız izleyici ölçümü ve jiroskop bias kestirimi ile MEKF.

Simülasyon; yörünge yayılımı, gerçek tutum üretimi, sensör ölçüm üretimi ve filtre sonuçlarının karşılaştırmalı görselleştirmesini içerir.

## Proje Akışı

`main.m` dosyası şu adımları yürütür:

1. `init.m` ile yörünge/sensör/filtre parametrelerini kurar.
2. `simulate_mode.m` ile coarse ve fine veri setlerini üretir.
3. Coarse başlangıcını `coarse_triad.m` ile hesaplar.
4. `coarse_mekf.m`, `eclipse_mekf.m` ve `mekf.m` ile filtreleri çalıştırır.
5. Quaternion ve Euler (RPY) tabanlı hata ölçümlerini ve grafiklerini üretir.

## Dosya Özeti

- `main.m`: Simülasyonun giriş noktası, filtre çağrıları, metrikler ve grafikler.
- `init.m`: Yörünge, sensör gürültüleri, bias ve kovaryans parametreleri.
- `simulate_mode.m`: Coarse/fine modlarında veri üretimi.
- `mekf.m`: Fine mode için MEKF (star tracker + gyro bias).
- `coarse_mekf.m`: Coarse mode için yön vektörü tabanlı MEKF.
- `eclipse_mekf.m`: Tutulma koşullarında gyro + manyetometre tabanlı MEKF.
- `coarse_triad.m`: Coarse başlangıç tutumunun TRIAD ile kestirimi.
- `simulate_*.m`: Sensör ölçüm üreticileri (sun sensor, magnetometer, gyro, star tracker).
- `propagate_orbit.m`, `attitude_truth.m`, `check_eclipse.m`: Yörünge, gerçek tutum ve tutulma modelleme yardımcıları.

## Gereksinimler

- MATLAB (depo MATLAB script/fonksiyonları ile çalışır)
- Ek toolbox bağımlılığı kod içinde zorunlu tanımlı değildir; temel MATLAB lineer cebir fonksiyonları kullanılır.

## Çalıştırma

Repo kök dizininde MATLAB komut penceresinde:

```matlab
main
```

## Üretilen Çıktılar

Çalışma sonunda:

- Coarse, eclipse ve fine mode için tutum karşılaştırma grafikleri (True vs Estimate),
- Yörünge grafiği,
- Sensör ölçüm grafikleri,
- Quaternion bileşen karşılaştırmaları,
- Komut penceresinde RMSE ve RMS/max hata özetleri üretilir.

## Lisans

Bu proje [MIT License](LICENSE) ile lisanslanmıştır.
