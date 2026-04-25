// Grup NO:12 Meryem Azra Gostak(240201061), Rana Huseynova (240201010)

#define _USE_MATH_DEFINES 
#define _CRT_SECURE_NO_WARNINGS

// Sabit RANSAC ayarlarý
#define MAX_ITERATIONS 1000
#define DISTANCE_THRESHOLD 0.02 // 2 cm tolerans (metre cinsinden)
#define MIN_INLIER_COUNT 8      // Proje isterine göre en az 8 nokta
#define MIN_LINE_SIZE 8         // RANSAC sonrasý minimum doðru büyüklüðü

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h> //SRAND icin
//#include <float.h> //INFINITY icin

// Struct tanýmý (sadece gerekli alanlar + çýktý olarak ranges eklendi)
typedef struct {
    double angle_min, angle_max, angle_increment;
    double range_min, range_max;
    double* ranges;   // filtrelenmiþ deðerler
    int     rangeCount; //gecerli olanlarin sayisi
    int toplamrange;//tum ranges degerlerinin sayisi
} lidarverileri;


typedef struct {
    double x;
    double y;
} kartezyenveri;

typedef struct {
    double A, B, C; // Doðru denklemi parametreleri
    int* inlierindisleri; // Bu doðruya ait noktalarýn orijinal dizideki indisleri
    int toplaminlier;
} dogru; // Doðruyu temsil eden yapý: Ax + By + C = 0

typedef struct {
    int max_deneme;         // Maksimum deneme sayýsý (Örn: 1000)
    double uzaklik_toleransi;  // Noktanýn doðruya uzaklýk toleransý (Örn: 0.05 metre)
    int min_gerekli_nokta;// Bir doðru kabul edilmesi için minimum nokta sayýsý (Projeye göre EN AZ 8)
} ransacayarlari;

typedef struct {
    double x;
    double y;
    int is_valid; // 1: Kesiþim var, 0: Paralel/Çakýþýk
} KesisimNoktasi;

// Ýmza ayný
lidarverileri verialma(FILE* dosya);  //verialma fonksiyonu tanimladik
kartezyenveri* kartezyene_cevir(const lidarverileri* l, int* nokta_sayisi); //kartezyenecevirme fonksiyonu tanimladik

dogru* dogrularin_ransacla_bulunmasi(const kartezyenveri* P, int nokta_sayisi, int* lineCount);// Ana fonksiyon imzasi
dogru dogru_hesaplama(const kartezyenveri P1, const kartezyenveri P2);// Ýki noktadan (P1, P2) doðru denklemini hesaplayan yardýmcý fonksiyon
double noktanin_dogruya_uzakligi(const kartezyenveri P, const dogru L);// Bir noktanýn (P) doðruya olan uzaklýðýný hesaplayan yardýmcý fonksiyon
KesisimNoktasi kesisim_hesapla(const dogru L1, const dogru L2); // Ýki doðrunun kesiþim noktasýný hesaplayan yardýmcý fonksiyon
double dogrular_arasi_aci(const dogru L1, const dogru L2); // Ýki doðru arasýndaki açýyý hesaplayan yardýmcý fonksiyon
double merkeze_mesafe_hesapla(const KesisimNoktasi P); //   Kesiþim noktasýnýn orijine (0,0) olan uzaklýðýný hesaplayan yardýmcý fonksiyon

// ---- SVG Çizim Prototipleri ----
void svg_ciz(
    const char* dosya_adi,
    const kartezyenveri* P, int N,
    const dogru* lines, int lineCount,
    const KesisimNoktasi* best_kesisim, // NULL olabilir
    double max_aci_deg,                 // 0 olabilir
    double mesafe_m                     // 0 olabilir
);

// ---------------- SVG Çizim Yardýmcýlarý ----------------
typedef struct { double minx, maxx, miny, maxy; } BBox;
static BBox bbox_hesapla(const kartezyenveri* P, int N);
static void world2svg(double x, double y, const BBox* b, int W, int H, int M, double* sx, double* sy);
static void dogruyu_kutu_ile_kes(const dogru* L, const BBox* b, double* x1, double* y1, double* x2, double* y2);
static void yaz_dogru_etiketi(FILE* f, int k, double sx1, double sy1, double sx2, double sy2);

int main() {

    srand(time(NULL)); // Rastgele sayý üreteciyi baþlat

    FILE* ScanData = fopen("C:\\Users\\Asus\\Desktop\\egitim\\UNIVERSITE\\2-1\\Programlama laboratuvari-1\\proje1 Lidar\\scan_data_NaN.toml", "r");
    if (ScanData == NULL) { printf("Dosya acilamadi!\n"); return 1; }
    printf("Dosya basariyla acildi.\n");

    lidarverileri lidar = verialma(ScanData);//altta yazilan fonksiyona atliyoruz
    fclose(ScanData);

    printf("\n--- Lidar Degerleri ---\n");
    printf("angle_min: %.7f\n", lidar.angle_min);
    printf("angle_max: %.7f\n", lidar.angle_max);
    printf("angle_increment: %.7f\n", lidar.angle_increment);
    printf("range_min: %.2f\n", lidar.range_min);
    printf("range_max: %.2f\n", lidar.range_max);
    printf("Filtrelenmis range sayisi: %d\n", lidar.rangeCount);
    printf("Toplam range sayisi: %d\n", lidar.toplamrange);
    //derece cinsinden yazma
    printf("angle_min[deg]: %.2f\n", lidar.angle_min * 57.29577951308232);
    printf("angle_max[deg]: %.2f\n", lidar.angle_max * 57.29577951308232);

    double	alinabilecek_max_deger_sayisi = (int)lround((lidar.angle_max - lidar.angle_min) / lidar.angle_increment) + 1;
    printf("Olmasi gereken deger sayisi: %5.2lf\n", alinabilecek_max_deger_sayisi);

    for (int i = 0; i < lidar.rangeCount; i++) {
        printf("%5.2lf ", lidar.ranges[i]);
        if ((i + 1) % 20 == 0) printf("\n");
    }
    printf("\n");

    // --- Kartezyene dönüþtür ---
    int nokta_sayisi = 0;
    kartezyenveri* noktalar = kartezyene_cevir(&lidar, &nokta_sayisi);

    // Tüm geçerli (NaN olmayan) noktalarý yazdýr
    printf("\n--- Gecerli noktalar (indis, aci[deg], r, x, y) ---\n");
    int yazilan = 0;
    for (int i = 0; i < nokta_sayisi; i++) {
        double r = lidar.ranges[i];
        double theta = lidar.angle_min + i * lidar.angle_increment;
        double deg = theta * 57.29577951308232;

        if (deg < lidar.angle_min * 57.29577951308232 || deg >alinabilecek_max_deger_sayisi)
        {
            continue;
        }
        // Sen geçersizleri 0.0 yapýyorsun; ayrýca kartezyene_cevir NaN yazýyor.
        // Ýki koþulu da güvenceye alalým:
        if (r > 0.0 && r >= lidar.range_min && r <= lidar.range_max
            && !isnan(noktalar[i].x) && !isnan(noktalar[i].y)) {

            printf("[%4d] %8.3f°  r=%7.3f -> x=%10.5f  y=%10.5f\n",
                i, deg, r, noktalar[i].x, noktalar[i].y);
            yazilan++;
        }
    }
    printf("Toplam gecerli nokta: %d\n", yazilan);

    int gecerli_nokta_sayisi = 0;
    for (int i = 0; i < nokta_sayisi; i++) {
        // Geçerli kabul edilen noktalarýn sadece NaN olmayanlar olmasý gerekir.
        // NaN kontrolü için isnan() fonksiyonunu kullanýyoruz.
        if (!isnan(noktalar[i].x) && !isnan(noktalar[i].y)) {
            gecerli_nokta_sayisi++;
        }
    }

    // RANSAC'a sadece geçerli noktalarý içeren yeni bir dizi gönderelim:
    kartezyenveri* RANSAC_noktalar = (kartezyenveri*)malloc(gecerli_nokta_sayisi * sizeof(kartezyenveri));
    int RANSAC_indis = 0;
    for (int i = 0; i < nokta_sayisi; i++) {
        if (!isnan(noktalar[i].x) && !isnan(noktalar[i].y)) {
            RANSAC_noktalar[RANSAC_indis++] = noktalar[i];
        }
    }

    printf("\n--- RANSAC Analizi ---\n");
    printf("RANSAC'a gonderilen gecerli nokta sayisi: %d\n", gecerli_nokta_sayisi);

    int lineCount = 0;
    // RANSAC ile doðrularý bul
    // Buradaki nokta_sayisi artýk gecerli_nokta_sayisi olmalýdýr!
    dogru* lines = dogrularin_ransacla_bulunmasi(RANSAC_noktalar, gecerli_nokta_sayisi, &lineCount);

    printf("Tespit edilen toplam dogru sayisi: %d\n", lineCount);
    printf("----------------------------------------\n");

    // Tespit edilen doðrularý ve inlier sayýlarýný listele
    for (int i = 0; i < lineCount; i++) {
        dogru L = lines[i];

        printf("\nDogru %d (Inlier: %d nokta):\n", i + 1, L.toplaminlier);

        // Doðru denklemi (Ax + By + C = 0)
        printf("  Denklem: %.4lf*x + %.4lf*y + %.4lf = 0\n", L.A, L.B, L.C);

        // Doðruya ait noktalarýn koordinatlarýný yazdýrma (Ýlk 5 nokta örnek olarak)
        printf("  Ornek Inlier Noktalar (ilk 5 - RANSAC dizisi indisleri):\n");
        int limit = (L.toplaminlier < 5) ? L.toplaminlier : 5;
        for (int j = 0; j < limit; j++) {
            // RANSAC_noktalar dizisinden koordinatlarý al
            int RANSAC_idx = L.inlierindisleri[j];
            printf("    (%d) x: %.3lf, y: %.3lf\n", RANSAC_idx, RANSAC_noktalar[RANSAC_idx].x, RANSAC_noktalar[RANSAC_idx].y);
        }
    }

    KesisimNoktasi best_kesisim = { 0 };
    double max_aci = 0.0;
    double final_mesafe = 0.0;
    int idx1 = -1, idx2 = -1;

    printf("\n--- Kesisim ve Aci Analizi ---\n");

    // Tüm doðru çiftlerini döngü ile kontrol et (i=d1, j=d2, ...)
    for (int i = 0; i < lineCount; i++) {
        for (int j = i + 1; j < lineCount; j++) {

            dogru L1 = lines[i];
            dogru L2 = lines[j];

            // 1. Kesiþim Noktasýný Bul
            KesisimNoktasi P = kesisim_hesapla(L1, L2);

            if (P.is_valid) {
                // 2. Doðrular Arasýndaki Açýyý Hesapla (Derece)
                double angle = dogrular_arasi_aci(L1, L2);

                printf("Dogru %d ve Dogru %d kesisiyor (%.2f, %.2f) - Aci: %.2f derece\n",
                    i + 1, j + 1, P.x, P.y, angle);

                // 3. Açý Kontrolü (60 derece ve üstü)
                if (angle >= 60.0) { // 60 derece, proje isterine göre 

                    // En yüksek açýlý kesiþimi bul (en uygun köþe)
                    if (angle > max_aci) {
                        max_aci = angle;
                        best_kesisim = P;
                        idx1 = i + 1;
                        idx2 = j + 1;
                    }
                }
            }
        }
    }

    // -------------------------------------------------------------------
    // ---- MERKEZ VE MESAFE HESABI (Ýster 5) ----
    // -------------------------------------------------------------------

    if (best_kesisim.is_valid) {
        // Robot merkezine olan mesafeyi hesapla
        final_mesafe = merkeze_mesafe_hesapla(best_kesisim); // Robot konumu (0, 0)

        printf("\n------------------------------------------------\n");
        printf("EN UYGUN HEDEF (Maksimum Aci: %.2f derece)\n", max_aci);
        printf("------------------------------------------------\n");
        printf("Kesen Doðrular: d%d ve d%d\n", idx1, idx2);
        printf("Kesiþim Noktasi (xc, yc): (%.4lf, %.4lf)\n", best_kesisim.x, best_kesisim.y);
        printf("Robot Merkezine Olan Mesafe: %.4lf metre\n", final_mesafe);
        printf("------------------------------------------------\n");
    }
    else {
        printf("\n60 derece ve üzeri açýyla kesiþen uygun bir doðru çifti bulunamadý.\n");
    }


    // iþ bitti
    free(noktalar);

    free(lidar.ranges);//tum islemler bittikten sonra yapcaz.aslinda kartezyene cevirdikten sonra yapilabilir.

    // --- SVG çýktý üret ---
    svg_ciz(
        "cikti.svg",
        RANSAC_noktalar, gecerli_nokta_sayisi,
        lines, lineCount,
        best_kesisim.is_valid ? &best_kesisim : NULL,
        max_aci,
        final_mesafe
    );
    printf("\nSVG olusturuldu: cikti.svg (tarayici ile acabilirsin)\n");

    return 0;
}

// ---- TOML oku + ranges'i basitçe ayýkla/filtrele ----
lidarverileri verialma(FILE* dosya) {
    lidarverileri lidar = { 0 };
    char line[1024]; //tampon(buffer) satir
    int inRanges = 0;                 // ranges = [ ... ] bloðunda mýyýz?
    double* gecerlidegerler = NULL; //bir double turunde dizi yaraticam fakat daha yaratmadim ama adresini kaydetmek icin pointer olusturdum.
    int n = 0;    // kac tane filtrelenmis yani gecerli deger var
    int y = 0; //toplam range degerlerinin sayisi

    while (fgets(line, sizeof(line), dosya)) {
        // 1) önce scaler parametreler
        if (inRanges == 0) {
            if (strstr(line, "angle_min"))            sscanf(line, "angle_min = %lf", &lidar.angle_min);
            else if (strstr(line, "angle_max"))       sscanf(line, "angle_max = %lf", &lidar.angle_max);
            else if (strstr(line, "angle_increment")) sscanf(line, "angle_increment = %lf", &lidar.angle_increment);
            else if (strstr(line, "range_min"))       sscanf(line, "range_min = %lf", &lidar.range_min);
            else if (strstr(line, "range_max"))       sscanf(line, "range_max = %lf", &lidar.range_max);
        }

        // 2) ranges baþlangýcý?
        if (inRanges == 0 && strstr(line, "ranges")) //inRanges ya da inRanges==NULL yazabiliriz
        {
            // bu satýrda dizi baþlayabilir
            if (strchr(line, '[')) //string icinde char arar.
                inRanges = 1;//range'in icinde oldugumuzu belirtir
        }

        // 3) ranges içindeysek token token topla
        if (inRanges == 1) {
            // virgül, köþeli parantez ve boþluklarý ayraç yap
            // Not: "ranges = [" olduðu satýrý da iþler; sonraki satýrlar da ayný mantýkla ilerler
            for (char* token = strtok(line, ",[] \t\r\n"); token; token = strtok(NULL, ",[] \t\r\n"))
            {
                // "ranges" ve "=" gibi baþlýk kýsýmlarýný atla
                if (strcmp(token, "ranges") == NULL || strcmp(token, "=") == NULL) continue;

                // sayýya çevir
                double alinanveri = atof(token); //ASCII to Float. Ranges kismini char olarak cektik ama ustunde kartezyene cevirme islemini yapmamiz icin double olmasi gerek.

                // filtre: range_min <= v <= range_max ve angle_min <= a <= angle_max

                if (alinanveri >= lidar.range_min && alinanveri <= lidar.range_max) {
                    // basitçe büyüt: (en basit hali — her eklemede realloc)
                    gecerlidegerler = (double*)realloc(gecerlidegerler, (n + 1) * sizeof(double));
                    gecerlidegerler[n++] = alinanveri;
                    y++;

                }
                else
                {
                    double alinanveri = 0;
                    gecerlidegerler = (double*)realloc(gecerlidegerler, (n + 1) * sizeof(double));
                    gecerlidegerler[n++] = 0.0;

                }
            }
            // Kapanýþ ayný satýrda mý? (']' varsa dizi bitti)
            if (strchr(line, ']')) inRanges = 0;
        }
    }

    lidar.ranges = gecerlidegerler;
    lidar.rangeCount = y;
    lidar.toplamrange = n;
    return lidar;
}

// ---- Polar(r,theta) -> Kartezyen(x,y) dönüþtürme ----
kartezyenveri* kartezyene_cevir(const lidarverileri* l, int* nokta_sayisi) {
    if (l == NULL || (*l).ranges == NULL || (*l).toplamrange <= 0) {
        if (nokta_sayisi) *nokta_sayisi = 0;
        return NULL;
    }

    int N = (*l).toplamrange;           // tüm ölçüm sayýsý (geçersizler dahil)
    kartezyenveri* P = (kartezyenveri*)malloc(N * sizeof(kartezyenveri));
    if (P == NULL) {
        if (nokta_sayisi) *nokta_sayisi = 0;
        return NULL;
    }

    for (int i = 0; i < N; i++) {
        double r = (*l).ranges[i];
        double theta = (*l).angle_min + i * (*l).angle_increment;

        // Geçerli mi? (Sen geçersizleri 0.0 olarak koyuyorsun)
        if (r > 0.0 && r >= (*l).range_min && r <= (*l).range_max) {
            P[i].x = r * cos(theta);
            P[i].y = r * sin(theta);
        }
        else {
            // Geçersizleri NaN iþaretleyelim
            P[i].x = NAN;
            P[i].y = NAN;
        }
    }

    if (nokta_sayisi) *nokta_sayisi = N;
    return P;
}

// Ýki noktadan (P1, P2) doðru denklemini hesaplayan yardýmcý fonksiyon
dogru dogru_hesaplama(const kartezyenveri P1, const kartezyenveri P2) {
    dogru L;
    L.A = P2.y - P1.y;  // A = y2 - y1
    L.B = P1.x - P2.x;  // B = x1 - x2 
    L.C = -(L.A * P1.x + L.B * P1.y);  // C = -(A*x1 + B*y1)

    L.inlierindisleri = NULL;
    L.toplaminlier = 0;

    return L;
}

// Bir noktanýn (P) doðruya olan dik uzaklýðýný hesaplayan yardýmcý fonksiyon
double noktanin_dogruya_uzakligi(const kartezyenveri P, const dogru L) {

    double pay = fabs(L.A * P.x + L.B * P.y + L.C); // Pay: |A*x0 + B*y0 + C|
    double payda = sqrt(L.A * L.A + L.B * L.B); // Payda: sqrt(A^2 + B^2)

    // Doðrunun sýfýr eðimli (A=B=0) olma durumunu kontrol et (genellikle RANSAC'ta olmaz, ama güvenlik için)
    if (payda < 1e-6) { // Çok küçük bir sayýya eþitse (sýfýr kabul et)
        return INFINITY; // Veya çok büyük bir hata deðeri
    }

    return pay / payda;
}

dogru* dogrularin_ransacla_bulunmasi(const kartezyenveri* P, int nokta_sayisi, int* lineCount) {
    if (nokta_sayisi < MIN_LINE_SIZE) {
        *lineCount = 0;
        return NULL;
    }

    // Doðrularý tutacak dinamik dizi
    dogru* lines = NULL;
    *lineCount = 0;

    // Noktalarýn kullanýlýp kullanýlmadýðýný takip etmek için bir bayrak dizisi
    // Baþlangýçta tüm noktalar KULLANILABÝLÝR (0)
    int* is_used = (int*)calloc(nokta_sayisi, sizeof(int));
    if (is_used == NULL) return NULL; // Bellek tahsis hatasý

    // Ýþlenmemiþ (kullanýlmamýþ) nokta sayýsýný takip et
    int remaining_points = nokta_sayisi;

    // Ana RANSAC döngüsü (Iterative RANSAC)
    // Kalan nokta sayýsý, bir doðru oluþturmak için yeterli olduðu sürece devam et
    while (remaining_points >= MIN_LINE_SIZE) {

        // Bu iterasyon için en iyi modeli tutacak deðiþkenler
        dogru bestLine = { 0, 0, 0, NULL, 0 };
        int bestInlierCount = 0;

        // -----------------------
        // RANSAC ÝÇ DÖNGÜSÜ BAÞLANGICI
        // -----------------------
        for (int i = 0; i < MAX_ITERATIONS; i++) {

            // 1. Rastgele Ýki Nokta Seçme (Benzersiz ve Kullanýlmamýþ olmalý)
            int idx1 = -1, idx2 = -1;
            int count_available = 0;

            // Sadece kullanýlmamýþ noktalarýn sayýsýný kontrol et
            for (int j = 0; j < nokta_sayisi; j++) {
                if (is_used[j] == 0) count_available++;
            }

            if (count_available < 2) break; // Ýki nokta seçilemiyorsa döngüyü kýr

            // Rastgele iki adet kullanýlmamýþ nokta seç
            do { idx1 = rand() % nokta_sayisi; } while (is_used[idx1] != 0);
            do { idx2 = rand() % nokta_sayisi; } while (is_used[idx2] != 0 || idx2 == idx1);

            kartezyenveri P1 = P[idx1];
            kartezyenveri P2 = P[idx2];

            // 2. Model Oluþturma (Doðru Denklemini Hesaplama)
            dogru currentLine = dogru_hesaplama(P1, P2);

            // 3. Uyum Kontrolü (Inlier'larý Bulma)
            int currentInlierCount = 0;
            int* currentInlierIndices = (int*)malloc(nokta_sayisi * sizeof(int));
            if (currentInlierIndices == NULL) continue; // Bellek hatasý

            for (int j = 0; j < nokta_sayisi; j++) {
                // Sadece kullanýlmamýþ noktalarý kontrol et
                if (is_used[j] == 0) {
                    double distance = noktanin_dogruya_uzakligi(P[j], currentLine);

                    // Mesafe tolerans içinde mi?
                    if (distance <= DISTANCE_THRESHOLD) {
                        currentInlierIndices[currentInlierCount++] = j;
                    }
                }
            }

            // 4. En Ýyi Modeli Kaydetme
            if (currentInlierCount > bestInlierCount) {
                // Önceki en iyi modelin inlier listesini serbest býrak
                if (bestLine.inlierindisleri != NULL) {
                    free(bestLine.inlierindisleri);
                }

                bestInlierCount = currentInlierCount;

                // Yeni en iyi modeli kaydet (Denklemi kopyala)
                bestLine.A = currentLine.A;
                bestLine.B = currentLine.B;
                bestLine.C = currentLine.C;

                // Inlier indislerini ve sayýsýný kaydet
                bestLine.inlierindisleri = currentInlierIndices;
                bestLine.toplaminlier = currentInlierCount;
            }
            else {
                // Mevcut model en iyi deðilse, inlier listesini serbest býrak
                free(currentInlierIndices);
            }
        }
        // -----------------------
        // RANSAC ÝÇ DÖNGÜSÜ BÝTÝÞÝ
        // -----------------------

        // 5. Model Kabulü ve Noktalarýn Çýkarýlmasý
        if (bestInlierCount >= MIN_LINE_SIZE) {
            // Yeni bir doðru bulundu! Doðruyu kaydet.
            (*lineCount)++;
            lines = (dogru*)realloc(lines, (*lineCount) * sizeof(dogru));
            if (lines == NULL) break; // Bellek hatasý

            // Doðruyu dinamik diziye ekle
            lines[(*lineCount) - 1] = bestLine;

            // **Çok Önemli:** Bu doðruya ait noktalarý "kullanýldý" olarak iþaretle.
            for (int k = 0; k < bestLine.toplaminlier; k++) {
                int index = bestLine.inlierindisleri[k];
                if (is_used[index] == 0) {
                    is_used[index] = 1; // Artýk kullanýlmýþ olarak iþaretle
                    remaining_points--;   // Kalan nokta sayýsýný azalt
                }
            }

        }
        else {
            // Yeterli inlier sayýsýna sahip bir doðru bulunamadý (veya kalan nokta azaldý)
            // Çýkmadan önce varsa bestLine.inlierIndices'i temizle
            if (bestLine.inlierindisleri != NULL) {
                free(bestLine.inlierindisleri);
            }
            break;
        }
    }

    // Bellek temizliði: Kullanýlan nokta bayrak dizisini serbest býrak
    free(is_used);

    // Tespit edilen doðrularý döndür
    return lines;
}

KesisimNoktasi kesisim_hesapla(const dogru L1, const dogru L2) {
    KesisimNoktasi P;
    P.is_valid = 0; // Baþlangýçta geçersiz

    // Diyagonal (Determinant) hesapla: D = A1*B2 - A2*B1
    double D = L1.A * L2.B - L2.A * L1.B;

    // Paralellik kontrolü (D yaklaþýk sýfýrsa, doðrular paraleldir)
    // 1e-6 gibi küçük bir tolerans deðeri kullanýlýr.
    if (fabs(D) < 1e-6) {
        // Doðrular paralel veya çakýþýk, kesiþim yok.
        return P;
    }

    // Kesiþim noktasýnýn koordinatlarýný hesapla
    // xc = (B1*C2 - B2*C1) / D
    // yc = (A2*C1 - A1*C2) / D
    P.x = (L1.B * L2.C - L2.B * L1.C) / D;
    P.y = (L2.A * L1.C - L1.A * L2.C) / D;

    P.is_valid = 1; // Kesiþim noktasý baþarýyla bulundu
    return P;
}

// Ýki doðru arasýndaki küçük açýyý (radyan cinsinden) hesaplar ve dereceye çevirir.
double dogrular_arasi_aci(const dogru L1, const dogru L2) {
    // Normal vektörler: n1 = (A1, B1) ve n2 = (A2, B2)
    double dot_product = L1.A * L2.A + L1.B * L2.B; // Nokta Çarpýmý
    double mag_n1 = sqrt(L1.A * L1.A + L1.B * L1.B); // Büyüklük 1
    double mag_n2 = sqrt(L2.A * L2.A + L2.B * L2.B); // Büyüklük 2

    // Güvenlik kontrolü
    if (mag_n1 < 1e-6 || mag_n2 < 1e-6) return 0.0;

    // Kosinüs hesapla: cos(theta) = |n1 . n2| / (|n1| |n2|)
    double cos_theta = fabs(dot_product) / (mag_n1 * mag_n2);

    // Kayan nokta hatalarýna karþý kontrol (cosinus -1 ile 1 arasýnda olmalý)
    if (cos_theta > 1.0) cos_theta = 1.0;
    if (cos_theta < 0.0) cos_theta = 0.0;

    // Açýyý radyan cinsinden bul ve dereceye çevir
    double angle_rad = acos(cos_theta);
    double angle_deg = angle_rad * 180.0 / M_PI;

    return angle_deg;
}

// Kesiþim noktasýnýn robotun konumuna (0, 0) olan mesafesini hesaplar (Ýster 5)
double merkeze_mesafe_hesapla(const KesisimNoktasi P) {

    return sqrt(P.x * P.x + P.y * P.y);  // Robot konumu (0, 0) olduðundan direkt Pisagor uygulanýr: d = sqrt(x^2 + y^2) [cite: 96]
}

void svg_ciz(const char* dosya_adi, const kartezyenveri* P, int N, const dogru* lines, int lineCount, const KesisimNoktasi* best_kesisim, double max_aci_deg, double mesafe_m)
{
    const int W = 1000, H = 800, M = 60; // geniþlik, yükseklik, kenar boþluðu
    FILE* f = fopen(dosya_adi, "w");
    if (!f) return;

    // Veri kutusu
    BBox box = bbox_hesapla(P, N);

    // Baþlýk
    fprintf(f, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
    fprintf(f, "<svg xmlns='http://www.w3.org/2000/svg' width='%d' height='%d'>\n", W, H);

    fprintf(f, "<svg xmlns='http://www.w3.org/2000/svg' width='%d' height='%d'>\n", W, H);
    fprintf(f, "<rect x='0' y='0' width='%d' height='%d' fill='white' stroke='none'/>\n", W, H);

    // Çerçeve
    fprintf(f, "<rect x='%d' y='%d' width='%d' height='%d' fill='none' stroke='black' stroke-width='1'/>\n",
        M, M, W - 2 * M, H - 2 * M);

    // Eksenler (robot merkezi (0,0))
    double sx0, sy0, sx1, sy1, sx2, sy2;
    world2svg(box.minx, 0, &box, W, H, M, &sx1, &sy1);
    world2svg(box.maxx, 0, &box, W, H, M, &sx2, &sy2);
    fprintf(f, "<line x1='%.2f' y1='%.2f' x2='%.2f' y2='%.2f' stroke='#888' stroke-width='1' />\n", sx1, sy1, sx2, sy2);
    world2svg(0, box.miny, &box, W, H, M, &sx1, &sy1);
    world2svg(0, box.maxy, &box, W, H, M, &sx2, &sy2);
    fprintf(f, "<line x1='%.2f' y1='%.2f' x2='%.2f' y2='%.2f' stroke='#888' stroke-width='1' />\n", sx1, sy1, sx2, sy2);

    // Noktalar
    for (int i = 0; i < N; ++i) {
        if (isnan(P[i].x) || isnan(P[i].y)) continue;
        double sx, sy; world2svg(P[i].x, P[i].y, &box, W, H, M, &sx, &sy);
        fprintf(f, "<circle cx='%.2f' cy='%.2f' r='2.2' fill='#1f77b4' />\n", sx, sy);
    }

    // Doðrular: inlier aralýðý kadar segment çiz
    for (int k = 0; k < lineCount; ++k) {
        const dogru* L = &lines[k];

        if (!L->inlierindisleri || L->toplaminlier < 2) {
            // Yedek: bbox ile kes (eski yöntem)
            double x1, y1, x2, y2; dogruyu_kutu_ile_kes(L, &box, &x1, &y1, &x2, &y2);
            world2svg(x1, y1, &box, W, H, M, &sx1, &sy1);
            world2svg(x2, y2, &box, W, H, M, &sx2, &sy2);
            fprintf(f, "<line x1='%.2f' y1='%.2f' x2='%.2f' y2='%.2f' stroke='#d62728' stroke-width='2' />\n",
                sx1, sy1, sx2, sy2);
            yaz_dogru_etiketi(f, k, sx1, sy1, sx2, sy2); // <-- EKLENDÝ

            continue;
        }

        // Yön vektörü: n = (A,B) normal ise, u = (-B, A)
        double ux = -L->B, uy = L->A;
        double un = sqrt(ux * ux + uy * uy);
        if (un < 1e-12) continue;
        ux /= un; uy /= un;

        double tmin = 0, tmax = 0;
        int init = 0;
        double x_min = 0, y_min = 0, x_max = 0, y_max = 0;

        for (int j = 0; j < L->toplaminlier; ++j) {
            int idx = L->inlierindisleri[j];
            double x = P[idx].x, y = P[idx].y;
            if (isnan(x) || isnan(y)) continue;
            double t = x * ux + y * uy;
            if (!init) { init = 1; tmin = tmax = t; x_min = x; y_min = y; x_max = x; y_max = y; }
            else {
                if (t < tmin) { tmin = t; x_min = x; y_min = y; }
                if (t > tmax) { tmax = t; x_max = x; y_max = y; }
            }
        }

        if (!init) continue;

        // (Ýstersen 2 cm kadar uzat)
        // x_min -= 0.02*ux; y_min -= 0.02*uy;
        // x_max += 0.02*ux; y_max += 0.02*uy;

        world2svg(x_min, y_min, &box, W, H, M, &sx1, &sy1);
        world2svg(x_max, y_max, &box, W, H, M, &sx2, &sy2);
        fprintf(f, "<line x1='%.2f' y1='%.2f' x2='%.2f' y2='%.2f' stroke='#d62728' stroke-width='2' />\n",
            sx1, sy1, sx2, sy2);
        yaz_dogru_etiketi(f, k, sx1, sy1, sx2, sy2); // <-- EKLENDÝ

    }


    // Robot merkezi
    world2svg(0, 0, &box, W, H, M, &sx0, &sy0);
    fprintf(f, "<circle cx='%.2f' cy='%.2f' r='5' fill='black'/>\n", sx0, sy0);
    fprintf(f, "<text x='%.2f' y='%.2f' font-size='14' font-family='Arial' fill='black'>Robot (0,0)</text>\n",
        sx0 + 8, sy0 - 8);

    // En iyi kesiþim + açý + mesafe
    if (best_kesisim && best_kesisim->is_valid) {
        double skx, sky; world2svg(best_kesisim->x, best_kesisim->y, &box, W, H, M, &skx, &sky);
        fprintf(f, "<circle cx='%.2f' cy='%.2f' r='5' fill='#2ca02c' />\n", skx, sky);
        fprintf(f, "<text x='%.2f' y='%.2f' font-size='14' font-family='Arial' fill='#2ca02c'>Kesisim (%.3f, %.3f)</text>\n",
            skx + 8, sky - 8, best_kesisim->x, best_kesisim->y);

        // Robot->Kesiþim mesafesi oku çiz
        fprintf(f, "<line x1='%.2f' y1='%.2f' x2='%.2f' y2='%.2f' stroke='#2ca02c' stroke-dasharray='5,5' stroke-width='1.5'/>\n",
            sx0, sy0, skx, sky);

        // Açý ve mesafe etiketi (üst bilgi)
        fprintf(f, "<text x='%d' y='%d' font-size='16' font-family='Arial' fill='black'>Aci: %.2f deg   Mesafe: %.3f m</text>\n",
            M + 6, M - 20 + 16, max_aci_deg, mesafe_m);

    }

    // Basit bir lejant
    int lx = W - M - 220, ly = M + 10;
    fprintf(f, "<rect x='%d' y='%d' width='210' height='90' fill='white' stroke='#ccc'/>\n", lx, ly);
    fprintf(f, "<circle cx='%d' cy='%d' r='3' fill='#1f77b4'/><text x='%d' y='%d' font-size='13' font-family='Arial'> Noktalar</text>\n",
        lx + 15, ly + 20, lx + 30, ly + 25);
    fprintf(f, "<line x1='%d' y1='%d' x2='%d' y2='%d' stroke='#d62728' stroke-width='2'/><text x='%d' y='%d' font-size='13' font-family='Arial'> Tespit edilen dogrular</text>\n",
        lx + 8, ly + 40, lx + 28, ly + 40, lx + 35, ly + 45);
    fprintf(f, "<circle cx='%d' cy='%d' r='5' fill='#2ca02c'/><text x='%d' y='%d' font-size='13' font-family='Arial'> Kesisim/ Hedef</text>\n",
        lx + 15, ly + 60, lx + 30, ly + 65);
    fprintf(f, "<circle cx='%d' cy='%d' r='5' fill='black'/><text x='%d' y='%d' font-size='13' font-family='Arial'> Robot (0,0)</text>\n",
        lx + 15, ly + 80, lx + 30, ly + 85);

    fprintf(f, "</svg>\n");
    fclose(f);
}

static BBox bbox_hesapla(const kartezyenveri* P, int N) {
    BBox b = { 0,0,0,0 };
    int ilk = 1;
    for (int i = 0; i < N; ++i) {
        if (isnan(P[i].x) || isnan(P[i].y)) continue;
        if (ilk) { b.minx = b.maxx = P[i].x; b.miny = b.maxy = P[i].y; ilk = 0; }
        else {
            if (P[i].x < b.minx) b.minx = P[i].x;
            if (P[i].x > b.maxx) b.maxx = P[i].x;
            if (P[i].y < b.miny) b.miny = P[i].y;
            if (P[i].y > b.maxy) b.maxy = P[i].y;
        }
    }
    // Hiç nokta yoksa, makul bir kutu ver
    if (ilk) { b.minx = -1; b.maxx = 1; b.miny = -1; b.maxy = 1; }
    // Biraz boþluk (padding)
    double dx = b.maxx - b.minx; if (dx <= 1e-9) dx = 1.0;
    double dy = b.maxy - b.miny; if (dy <= 1e-9) dy = 1.0;
    double pad = 0.10; // %10
    b.minx -= dx * pad; b.maxx += dx * pad;
    b.miny -= dy * pad; b.maxy += dy * pad;
    // (0,0) robot merkezini kadraja dahil et
    if (0 < b.minx) b.minx = -0.1;
    if (0 > b.maxx) b.maxx = 0.1;
    if (0 < b.miny) b.miny = -0.1;
    if (0 > b.maxy) b.maxy = 0.1;
    return b;
}

static void world2svg(double x, double y, const BBox* b, int W, int H, int M, double* sx, double* sy)
{
    // SVG sol-üst (0,0), y aþaðý artar. Dünyayý geniþliðe/yeleðe orantýla.
    double sx_scale = (W - 2.0 * M) / (b->maxx - b->minx);
    double sy_scale = (H - 2.0 * M) / (b->maxy - b->miny);
    double s = sx_scale < sy_scale ? sx_scale : sy_scale; // eþ ölçek
    double cx = M + (x - b->minx) * s;
    double cy = M + (b->maxy - y) * s; // y eksenini ters çevir
    *sx = cx; *sy = cy;
}

static void dogruyu_kutu_ile_kes(const dogru* L, const BBox* b, double* x1, double* y1, double* x2, double* y2)
{
    // Doðru: A x + B y + C = 0
    // Görünürlük için kutu sýnýrlarýyla kesiþtirip 2 nokta üretelim.
    // Kutu kenarlarý: x = minx/maxx, y = miny/maxy
    // Olasý 4 kesiþimden ekran içinde kalan 2'sini seç.
    double xs[4], ys[4]; int c = 0;

    // x = minx => y = (-A*x - C)/B  (B ~ 0 ise atla)
    if (fabs(L->B) > 1e-12) {
        double y = (-(L->A * b->minx) - L->C) / L->B;
        if (y >= b->miny && y <= b->maxy) { xs[c] = b->minx; ys[c] = y; c++; }
    }
    // x = maxx
    if (fabs(L->B) > 1e-12) {
        double y = (-(L->A * b->maxx) - L->C) / L->B;
        if (y >= b->miny && y <= b->maxy) { xs[c] = b->maxx; ys[c] = y; c++; }
    }
    // y = miny => x = (-B*y - C)/A  (A ~ 0 ise atla)
    if (fabs(L->A) > 1e-12) {
        double x = (-(L->B * b->miny) - L->C) / L->A;
        if (x >= b->minx && x <= b->maxx) { xs[c] = x; ys[c] = b->miny; c++; }
    }
    // y = maxy
    if (fabs(L->A) > 1e-12) {
        double x = (-(L->B * b->maxy) - L->C) / L->A;
        if (x >= b->minx && x <= b->maxx) { xs[c] = x; ys[c] = b->maxy; c++; }
    }

    if (c >= 2) {
        *x1 = xs[0]; *y1 = ys[0];
        *x2 = xs[1]; *y2 = ys[1];
    }
    else {
        // Kutu içinde yeterli kesiþim yoksa, merkezden kýsa bir segment çiz
        *x1 = b->minx; *y1 = (-(L->A * (*x1)) - L->C) / (fabs(L->B) < 1e-12 ? 1.0 : L->B);
        *x2 = b->maxx; *y2 = (-(L->A * (*x2)) - L->C) / (fabs(L->B) < 1e-12 ? 1.0 : L->B);
    }
}

// Ekran koordinatlarýndaki (sx1,sy1)-(sx2,sy2) çizgisine "dK" etiketi yazar
static void yaz_dogru_etiketi(FILE* f, int k, double sx1, double sy1, double sx2, double sy2) {
    double dx = sx2 - sx1, dy = sy2 - sy1;
    double L = sqrt(dx * dx + dy * dy);
    if (L < 1e-6) return;

    // Orta nokta + çizgiden 10 px yana (normal doðrultusunda) kaydýr
    double mx = 0.5 * (sx1 + sx2), my = 0.5 * (sy1 + sy2);
    double nx = -dy / L, ny = dx / L;   // ekrana göre normal
    mx += 10.0 * nx;
    my += 10.0 * ny;

    // Yazýyý çizgiye paralel döndür
    double angle_deg = atan2(dy, dx) * 180.0 / M_PI;

    char etiket[16];
    sprintf(etiket, "d%d", k + 1);

    // Beyaz konturlu kýrmýzý yazý (arkasý okunaklý olsun diye)
    fprintf(f,
        "<text x='%.2f' y='%.2f' font-size='13' font-family='Arial' "
        "fill='#d62728' transform='rotate(%.2f, %.2f, %.2f)' "
        "style='paint-order:stroke; stroke:white; stroke-width:3px'>%s</text>\n",
        mx, my, angle_deg, mx, my, etiket);

    // Ýnce bir üst yazý (keskinlik için)
    fprintf(f,
        "<text x='%.2f' y='%.2f' font-size='13' font-family='Arial' "
        "fill='#d62728' transform='rotate(%.2f, %.2f, %.2f)'>%s</text>\n",
        mx, my, angle_deg, mx, my, etiket);
}