# hse21_hw1
Ссылка на Google Colab, содержащий код для анализа контигов и скаффолдов: https://colab.research.google.com/drive/19yNnDD2RGLE6bUPGIHZmgD30wg8aOJL0?usp=sharing.

Создание символических ссылок на файлы
```
ln -s /usr/share/data-minor-bioinf/assembly/oil_R1.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oil_R2.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R1_001.fastq
ln -s /usr/share/data-minor-bioinf/assembly/oilMP_S4_L001_R2_001.fastq
```
Выбор случайных чтений
```
seqtk sample -s928 oil_R1.fastq 5000000 > oil_R1_sample.fastq
seqtk sample -s928 oil_R2.fastq 5000000 > oil_R2_sample.fastq
seqtk sample -s928 oilMP_S4_L001_R1_001.fastq 1500000 > oilMP_S4_L001_R1_001_sample.fastq
seqtk sample -s928 oilMP_S4_L001_R2_001.fastq 1500000 > oilMP_S4_L001_R2_001_sample.fastq
```
Удаление символических ссылок
```
rm oil_R1.fastq
rm oil_R2.fastq
rm oilMP_S4_L001_R1_001.fastq
rm oilMP_S4_L001_R2_001.fastq
```
Оценка качества исходных чтений (выдача производится в отдельную папку fastqc)
```
mkdir fastqc
ls *.fastq | xargs -P 4 -tI{} fastqc -o fastqc {}
```
Оъединение файлов при помощи multiqc
```
mkdir multiqc
multiqc -o multiqc fastqc
```
Подрезание чтений по качеству и удаление праймеров
```
platanus_trim oil_R1_sample.fastq oil_R2_sample.fastq
platanus_internal_trim oilMP_S4_L001_R1_001_sample.fastq oilMP_S4_L001_R2_001_sample.fastq
```
Удаление исходных неподрезанных файлов
```
rm oilMP_S4_L001_R1_001_sample.fastq
rm oil_R2_sample.fastq
rm oilMP_S4_L001_R2_001_sample.fastq
rm oil_R1_sample.fastq
```
Перемещение подрезанных файлов в отдельную папку
```
mkdir trimmed_fastq
mv -v *trimmed trimmed_fastq/
```
Оценка качества подрезанных чтений
```
mkdir trimmed_fastqc
ls trimmed_fastq/* | xargs -P 4 -tI{} fastqc -o trimmed_fastqc {}
```
Объединение подрезанных файлов
```
mkdir trimmed_multiqc
multiqc -o trimmed_multiqc trimmed_fastqc
```
Сборка контигов при помощи platanus assemble
```
time platanus assemble -o Poil -t 1 -m 25 -f trimmed_fastq/oil_R1.fastq.trimmed trimmed_fastq/oil_R2.fastq.trimmed 2> assemble.log
```
Сборка скаффолдов при помощи platanus scaffold
```
time platanus scaffold -o Poil -t 1 -c Poil_contig.fa -IPl trimmed_fastq/oil_R1_sample.fastq.trimmed trimmed_fastq/oil_R2_sample.fastq.trimmed -OP2 trimmed_fastq/oilMP_S4_L001_R1_001_sample.fastq.int_trimmed trimmed_fastq/oilMP_S4_L001_R2_001_sample.fastq.int_trimmed 2> scaffold.log
```
Поиск самого длинного скаффолда
```
head Poil_scaffold.fa
```
Данный скаффолд называется >scaffold1_len3838373_cov231

Сохраняем отдельно данный файл
```
echo scaffold1_len3838373_cov231 > _tmp.txt
seqtk subseq Poil_scaffold.fa _tmp.txt > scaffold1_len3838373_cov231.fasta
```
Подсчет числа гэпов
```
tail -n +2 scaffold1_len3838373_cov231.fasta | perl -pe 's/[ACTG]+/\n/gi;' | grep -e 'N' | wc -l
```
Получилось 65 гэпов

Подсчет длины гэпов
```
tail -n +2 scaffold1_len3838373_cov231.fasta | perl -pe 's/[ACTG]+//gi;' | wc
```
Длина гэпов = 6347

Уменьшение количества гэпов при помощи gap_close
```
time platanus gap_close -o Poil -t 1 -c Poil_scaffold.fa -IPl trimmed_fastq/oil_R1_sample.fastq.trimmed trimmed_fastq/oil_R2_sample.fastq.trimmed -OP2 trimmed_fastq/oilMP_S4_L001_R1_001_sample.fastq.int_trimmed trimmed_fastq/oilMP_S4_L001_R2_001_sample.fastq.int_trimmed 2> gapclose.log
```
Снова находим самый длинный скаффолд в подрезанных чтениях
```
head Poil_gapClosed.fa
```
Это scaffold1_cov231.
Отдельно сохраняем данный скаффолд
```
echo scaffold1_cov231 > _tmp.txt
seqtk subseq Poil_gapClosed.fa _tmp.txt > oil_genome.fna
```
Считаем измененное число гэпов
```
tail -n +2 oil_genome.fna | perl -pe 's/[ACGT\n\r]+/\n/gi;' | grep -e 'N' | wc -l
```
Осталось 9 гэпов

Считаем длину гэпов
```
tail -n +2 oil_genome.fna | perl -pe 's/[ACGT\n\r]+//gi;' | wc
```
Длина оставшихся гэпов = 1674

Скриншоты multiqc до подрезания чтений
<img width="717" alt="general_statistics_multiqc" src="https://user-images.githubusercontent.com/91226664/139124464-89b3bb08-01cc-474c-b44a-af666c52618a.png">
<img width="698" alt="mean_quality_scores_multiqc" src="https://user-images.githubusercontent.com/91226664/139124922-df41c821-bd5a-466e-8dad-1d865608aaf4.png">
<img width="698" alt="per_sequence_quality_scores_multiqc" src="https://user-images.githubusercontent.com/91226664/139125217-d5784928-e4a5-4265-b38b-9f4b418b7450.png">
<img width="689" alt="adapter_content_multiqc" src="https://user-images.githubusercontent.com/91226664/139125416-382ac40d-b30f-49a4-8239-a385fcd67325.png">
Скриншоты multiqc после подрезания чтений
<img width="721" alt="general_statistics_multiqc_trimmed" src="https://user-images.githubusercontent.com/91226664/139125784-aa346d7b-2e68-4998-a629-18e8ddded021.png">
<img width="693" alt="mean_quality_scores_multiqc_trimmed" src="https://user-images.githubusercontent.com/91226664/139126321-0c238a95-c14e-463f-9d0d-e0b4c99e0b5c.png">
<img width="696" alt="per_sequence_quality_scores_multiqc_trimmed" src="https://user-images.githubusercontent.com/91226664/139126332-3b375547-1567-4a52-9e0c-db556c786cee.png">
<img width="695" alt="adapter_content_multiqc_trimmed" src="https://user-images.githubusercontent.com/91226664/139126343-2c7ac375-e8eb-4b83-bc2c-a5eeeef3ae44.png">

