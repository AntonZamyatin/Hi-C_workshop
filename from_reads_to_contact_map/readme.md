Пошаговая инструкция для составления пайплайна на основе [Arima pipeline](https://github.com/ArimaGenomics/mapping_pipeline/blob/master/Arima_Mapping_UserGuide_A160156_v02.pdf)

[Vertebrate Genome Project (VGP) версия пайплайна](https://github.com/VGP/vgp-assembly/blob/master/pipeline/salsa/arima_mapping_pipeline.sh)

## Данные
На вход пайплайна поступают **два .fastq файла с прямыми и обратными ридами секвенирования библиотеки Hi-C**.
Также необходим **референсный .fasta** файл с референсными последовательностями или геномом, на который будет осуществлятся мэппинг.

для самостоятельной работы скачайте референс генома человека hg38 в fasta формате:
https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/
и риды Hi-C:

[forward](https://s3.amazonaws.com/4dn-dcic-public/hic-data-analysis-bootcamp/input_R1.fastq.gz)

[reverse](https://s3.amazonaws.com/4dn-dcic-public/hic-data-analysis-bootcamp/input_R2.fastq.gz)

Эти данные взяты с репозитория Hi-C bootcamp с хорошими 
описаниями пайплайнов в виде презентаций c воркшопов. 
Ссылка: https://github.com/hms-dbmi/hic-data-analysis-bootcamp

## Настройка окружения
Рекомендуется  провести предварительную подготовку окружения до того, как запускать пайплайн.
Одним из хороших способов подготовки окружения является использование [miniconda3](https://docs.conda.io/en/latest/miniconda.html) менеджера пакетов иокружений.
После установки Conda рекомендуется установить более современный пакетный менеджер mamba. создать чистое окружение и активировать его. 
```bash
conda install mamba
conda create -n your_environment_name
conda activate your_environment_name
```
Установка инструментов необходимых в пайплайне:
1. [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
```bash
mamba install -c bioconda fastqc
```
2. [BWA](https://github.com/lh3/bwa)
```bash
mamba install -c bioconda bwa
```
3. [Samtools](https://www.htslib.org/)
```bash
mamba install -c bioconda samtools
```
4. Скрипты Arima  пайплайна на Perl:
   * filter_5end.pl
   * two_read_bam_combiner.pl
   * get_stats.pl
5. Скомпилированный .jar файл [Picard](https://github.com/broadinstitute/picard) его обязательно запускать на Java выше 17ой версии. Скачать скомпилированный picard.jar с последнего релиза можно тут: 
https://github.com/broadinstitute/picard/releases/tag/3.0.0

6. create_chrom_sizes.py Скрипт для создания файла с размерами хромосом референсаю 
7. [Cooler](https://github.com/open2c/cooler)
```bash
mamba install -c bioconda cooler
```

## 1. Контроль качества секвенирования и фильтрация ридов если необходимо
```bash
fastqc -t ${CPU} reads_R1.fastq reads_R2.fastq -o outdir
```

## 2. Индексация референса
Используем bwa для получения файлов индекса: .amb .ann .bwt .pac .sa
```bash
bwa index reference.fasta
```

## 3. Выравнивание прямых и обратных ридов на референс
Создаем директорию для первичного выравнивания:
```bash
mkdir RAW_ALN
```
Выравниваем по отдельности прямые и обратные ряды.
```
bwa mem -t ${CPU} reference.fasta forward_reads.fastq |\
samtools view -Sb > RAW_ALN/R1_aln.bam

bwa mem -t ${CPU} reference.fasta reverse_reads.fastq |\
samtools view -Sb > RAW_ALN/R2_aln.bam
```
**bwa mem** выполняет непосредственное выравнивание и сохраняет его в текстовом SAM формате.
**samtools view -Sb** переводит SAM формат в двоичный BAM

## 4. 5' фильтрация
создаем директорию для фильтрованных выравниваний
```bash
mkdir FILT_ALN
```
используем **filter_5end.pl** из Arima пайплайна для фильтрации каждого файла выравнивания по отдельности.
```bash
samtools view -@ ${CPU} -h RAW_ALN/R1_aln.bam |\
perl filter_5end.pl |\
samtools view -@ ${CPU} -Sb  > FILT_ALN/R1_aln_filt.bam

samtools view -@ ${CPU} -h RAW_ALN/R2_aln.bam |\
perl filter_5end.pl |\
samtools view -@ ${CPU} -Sb  > FILT_ALN/R2_aln_filt.bam
```
здесь BAM переводитсяв текстовый формат, обрабатывается скриптом фильтрации и запаковывается обратно в двоичный формат

## 5. Соединяем выравнивания
создаем директорию для комбинированных выравниваний
```bash
mkdir COMB_ALN
```
используем скрипт **two_read_bam_combiner.pl** из Arima  пайплайна
```bash
perl two_read_bam_combiner.pl \
	FILT_ALN/R1_aln_filt.bam \
	FILT_ALN/R2_aln_filt.bam |\
samtools view -@ ${CPU} -Sb > COMB_ALN/aln.bam
```

## 6. Дедупликация
создаем директорию для дедуплицированных выравниваний и директорию для временных файлов
```bash
mkdir DEDUP_ALN
mkdir DEDUP_ALN/tmp
```
сортируем выравнивания по координатам:
```bash
samtools sort \
	-@ ${CPU} \
	-T DEDUP_ALN/sort_aln.bam.tmp \
	-m2G \
	-O bam \
	-o DEDUP_ALN/sort_aln.bam \
	COMB_ALN/aln.bam
```
дедупликация с помощью Picard:
```bash
java -jar -Xmx6g -Djava.io.tmpdir=DEDUP_ALN/tmp picard.jar MarkDuplicates \
	-REMOVE_DUPLICATES true \
	-I DEDUP_ALN/sort_aln.bam \
	-O DEDUP_ALN/dedup_aln.bam \
	-M DEDUP_ALN/dedup_aln.metrics.txt \
	-ASSUME_SORT_ORDER coordinate \
	-MAX_FILE_HANDLES_FOR_READ_ENDS_MAP 1024 
```
пересортировываем дедуплицированные выравнивания в порядке имен:
```bash
samtools sort \
	-@ ${CPU} \
	-n \
	-T DEDUP_ALN/dedup_aln_sortbyname.bam.tmp \
	-m2000m \
	-O bam \
	-o DEDUP_ALN/dedup_aln_sortbyname.bam \
	DEDUP_ALN/dedup_aln.bam
```

## 7. Получение статистики
используем скрипт **get_stats** из Arima пайплайна
```bash
perl get_stats.pl DEDUP_ALN/dedup_aln_sortbyname.bam > \
	DEDUP_ALN/dedups_bam.stats
```
на этом закончена часть с выравниваниями и начинается построение карты Hi-C контактов при помощи cooler

## 8. Подготовка файла с размерами хромосом
Создаем директорию для работы с Cooler
```bash
mkdir COOLER
```
Необходимо получить текстовый файл с размерами хромосом в формате:
	chr1    248956422
	chr2    242193529
	chr3    198295559
	chr4    190214555
	chr5    181538259
	...
Такой файл можно составить самим или использовать скрипт **create_chrom_sizes.py**. Для его работы неоходимо установить biopython
```bash
mamba install biopython
python create_chrom_sizes.py reference.fasta out_directory
```

## 9. Перевод файла выравнивания в текстовый формат таблицы контактов
используем **bam2pairs**,данный инструмент входит в состав пакета cooler
```bash
label='your_name_for_output'
bam2pairs -c COOLER/chrom.sizes DEDUP/dedup_aln_sortbyname.bam COOLER/${label}
```

## 10. Биннинг с помощью cooler
Используем таблицу контактов для создания карты контактов с одним разрешение, с указанием размера бина (длина в парах нуклеотидов каждого пикселя карты контактов).
Для создания .cool файла с единственным разрешением используем команду **cload pairix**
```bash
cooler cload pairix \
	-p ${CPU} \
	COOLER/chrom.sizes:100000 \
	COOLER/${label}.bsorted.pairs.gz \
	COOLER/${label}_100k.cool
```

## 11. Балансировка и создание дополнительных слоев разрешений матрицы
Последним этапом получения матрицы контактов будет балансировка и создание более низких разрешений для удобства визуализации.
для этого будем использовать **cooler zoomify**
```bash
cooler zoomify \
	-n ${CPU} \
	-r 100000N, \
	--balance \
	--balance-args "--nproc ${CPU}" \
	-o COOLER/${label}_multires.mcool \
	COOLER/${label}_100k.cool
```
\-r - опция какие разрешения применить 100000N оббозначает использовать разрешения 100000, 200000, 500000, 1000000...
есть также опция 4DN для карт, изначально полученных в разрешении 1000 (самое высокое)
4DN (4D nucleome) - организация, которая ввела удобный стандарт разрешений для Hi-C карт: 1000, 2000, 5000, 10000, 25000 ...

На этом этапе мы получили конечную карту контактов, с которой можем работать дальше. 
Чтобы её визуализировать можно использовать [HiGLass local server](https://github.com/higlass/higlass), который можно поставить через докер и утилиту **higlass-manager** (смотрите документацию).
Также можно использовать API самого Cooler, ноутбук с описанием API находится здесь: https://github.com/open2c/cooler-binder/blob/master/cooler_api.ipynb