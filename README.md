# Домашняя работа №4
### Пайплайн bioinfo_joint для обработки последовательностей РНК и ДНК 


Пайплайн содержит биоинформатические утилиты, которые позволяют транслировать, транскрибировать, обращать и обратно транскрибировать последовательности ДНК, а также фильтровать наборы последовательностей ДНК по GC-составу, качеству прочтений и длине.

## Возможности проекта

1. Основной скрипт bioinfo_joint содержит утрилиту для проведения операций над ДНК и РНК. Утилита принимает на вход строки, состоящие из последовательности аминокислот, содержащие в последней позиции команду-операцию, которую необходимо применить к последовательностям нуклеиновых кислот. Операции включают `translate` (трансляция), `transcribe` (транскрипция), `complement` (комплементарная последовательность), `reverse_complement` (обратная комплементарная последовательность). В ответ на ввод последовательностей, не соотвествующих ДНК или РНК, выводится 0.

2. Также основной скрипт bioinfo_joint содержит утрилиту для отбора последовательностей из словаря по длине, CG-составу и качеству прочтений. На вход следует подавать словарь с последовательностями прочтений fastq и сроку с качеством эти прочтений, задать интервал GC состава (в процентах), интервал длины и пороговое значение качества по шкале phred33. На выходе вторая утилита возвращать аналогичный словарь, состоящий только из тех сиквенсов, которые удовлетвопмли всем условиям:
</br></br>


## Кому может быть полезен этот проект

Проект может быть полезен молекулярным биологам для планирования кспериментов in silico, в частности, при работе с последовательносяти ДНК и РНК для создания праймеров и предсказания транскриптов. Также пакет может быть использован для фильтрации данных, полученных с помощью технологий полногеномного секвинирования.

## Кому может быть полезен этот проект

Чтобы начать работать с проектом необходимо импортировать главный модуль bioinfo_joint и использовать его целиком или отдельные функции из модуля для обработки своих данных. Кроме основных функций `run_dna_rna_tools` и `filter_fastq` следует отметить еще несколько полезных функций проекта, которые можно использовать, как самостоятельные, :
`translate`
`reverse`
`complement`
`reverse_complement`
`qc_calc`
`gc_cont`.

## Помощь и обратная связь по проекту
Чтобы обратиться за помощью по проекту, можно связаться напрямую с автором проекта по электронной почте elena.n.kozhevnikova@gmail.com или можно написать в репозитории в Issue Tracker и оставить свой вопрос. На него с радостью ответит очень милая девушка :information_desk_person:


## Авторство и вклад в проект
Проект создан Elena-Kozhevnikova (https://github.com/Elena-Kozhevnikova) на основе учебной практики в Институте Биоинформатики :heartpulse: с помощью и поддержкой команды института и особенно Nikita Vaulin (https://github.com/nvaulin).
