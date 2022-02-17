# hse_hw1_meth

Ссылка на код: https://colab.research.google.com/drive/1f8q550GayvWvzzW3RV7DyIf0-FazO5Re?usp=sharing
### FQ-анализ (выполнен для образца 8 Cell)

![](https://github.com/kolbunovaa/images/blob/main/2022-02-17_00-57-34.png)

![](https://github.com/kolbunovaa/images/blob/main/2022-02-17_00-58-07.png)

### FQ-анализ для RNA-seq образца из одного из прошлых дз

![](https://github.com/kolbunovaa/images/blob/main/2022-02-17_13-33-25.png)

![](https://github.com/kolbunovaa/images/blob/main/2022-02-17_13-33-50.png)

Для секвенирования ДНК наблюдается примерно, с небольшими отклонениями, одинаковое содержание каждого нуклеотида (правда в начале ридов в данном образце есть проблемы), тогда как для бисульфитного секвенирования наблюдается разброс в содержании каждого из нуклеотидов. Так,  уровень содержания A, G приблизительно одинаков, T имеет повышенное содержание, а C - пониженное. Так как FQ-анализ проводился для 8 Cell (для которых уровень метилирования не очень высокий), то картина в целом выглядит логично. Неметилированные цитозины при данном секвенировании превращаются в урацилы (которым комплементарны тимины), метилированные - остаются неизменными. Так как уровень метилирования небольшой, то и уровень неизменных цитозинов будет низким, а количество тиминов возрастет.


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

Первое, что бросается в глаза, - это то, что CHG- и  CHH-метилирования имеют очень низкий процент во всех образцах (для дальних позициях в ридах они практически везде на уровне 0%), основным является CpG-метилирование. Самый высокий уровень метилирования в образце Epiblast, почти 80%; самый низкий - в ICM (примерно 25%); в образце 8 Cell - находится примерно на уровне 45%.

e)
Распределение метилирования

![](https://github.com/kolbunovaa/images/blob/main/8cell.png)

![](https://github.com/kolbunovaa/images/blob/main/epi.png)

![](https://github.com/kolbunovaa/images/blob/main/icm.png)

Для 8 Cell: большой пик при 0%, однако есть пара существенных пиков при 50% и 100%, и в целом они друг друга компенсируют, также есть несколько пиков в середине графика, которые образуют кривую, очень похожую на распределение Гауса; в целом можно сказать, что метилирование находится на среднем уровне.

Самый высокий пик среди всех образцов имеет ICM при значении 0% метилирования, остальные небольшие пики в основном находятся в области 0-50%, а пик на 100% имеет довольно маленькую частоту (примерно 0,05), это свидетельствует о низком метилировании цитозинов.

Для Epiblast основная масса пиков находится в высокопроцентной области, пик на 100% почти достигает частоты 0,5; что говорит о высоком уровне метилирования.

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

