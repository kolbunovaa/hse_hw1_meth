# hse_hw1_meth

Ссылка на код: https://colab.research.google.com/drive/1f8q550GayvWvzzW3RV7DyIf0-FazO5Re?usp=sharing
### FQ-анализ (выполнен для 8 Cell образца)

![](https://github.com/kolbunovaa/images/blob/main/2022-02-17_00-57-34.png)

![](https://github.com/kolbunovaa/images/blob/main/2022-02-17_00-58-07.png)

a)
Сводная таблица по количествам ридов 

|    Region    |     8 Cell     |    Epiblast   |    ICM     | 
| :---         |     :---:      |     :---:     |    ---:    | 
| 11347700-11367700   | 1090     | 2328    | 1456       |              
| 40185800-40195800   | 464      | 1062    |  630       | 

b)
Уровень дупликации для каждого из образцов
|     | 8 Cell   | Epiblast    |  ICM  |
| :--- | :---: | :---: | ---:   |
|Duplication, % | 18.31  | 2.92  | 9.08  |

*bash-скрипт: ! ls *pe.bam | xargs -P 4 -tI{} deduplicate_bismark  --bam  --paired  -o s_{} {}

d) M-bias plot

## 8 Cell

![](https://github.com/kolbunovaa/images/blob/main/8celBismark%20M-bias%20Read%201.png)

![](https://github.com/kolbunovaa/images/blob/main/8cellBismark%20M-bias%20Read%202%20(1).png)

## Epiblast

![](https://github.com/kolbunovaa/images/blob/main/epi_Bismark%20M-bias%20Read%201.png)

![](https://github.com/kolbunovaa/images/blob/main/epi_Bismark%20M-bias%20Read%202%20(1).png)

## ICM

![](https://github.com/kolbunovaa/images/blob/main/icm_Bismark%20M-bias%20Read%201.png)

![](https://github.com/kolbunovaa/images/blob/main/icm_Bismark%20M-bias%20Read%202.png)

e)
Распределение метилирования

![](https://github.com/kolbunovaa/images/blob/main/8cell.png)

![](https://github.com/kolbunovaa/images/blob/main/epi.png)

![](https://github.com/kolbunovaa/images/blob/main/icm.png)

Код:
```
import pandas as pd

first = pd.read_csv("s_SRR5836473_1_bismark_bt2_pe.deduplicated.bedGraph", delimiter='\t', skiprows=1, header=None)
first.head()

second = pd.read_csv("s_SRR3824222_1_bismark_bt2_pe.deduplicated.bedGraph", delimiter='\t', skiprows=1, header=None)
second.head()

third = pd.read_csv("s_SRR5836475_1_bismark_bt2_pe.deduplicated.bedGraph", delimiter='\t', skiprows=1, header=None)
third.head()
```
```
import matplotlib.pyplot as plt

plt.figure(figsize=[11, 8])
plt.title("Распределение метилирования в 8 Cell", fontsize=19)
plt.hist(first[3], bins=100, density=True)
plt.show()

plt.figure(figsize=[11, 8])
plt.title("Распределение метилирования в Epiblast", fontsize=19)
plt.hist(second[3], bins=100, density=True)
plt.show()

plt.figure(figsize=[11, 8])
plt.title("Распределение метилирования в ICM", fontsize=19)
plt.hist(third[3], bins=100, density=True)
plt.show()
```

f)
Визуализация метилирования и покрытия

![](https://github.com/kolbunovaa/images/blob/main/image_cov2.png)


![](https://github.com/kolbunovaa/images/blob/main/image_cov%20(1).png)

