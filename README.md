# Readme

## Что сделано

Реализация эллиптической кривой в форме Гессе и арифметических операций над точкой, с использованием библиотеки gmp.<br/>

- Операция сложения реализована алгоритмом “add-2009-bkl”.<br/>
- Операция удвоения реализована алгоритмом “dbl”.<br/>
- Для вычисления кратной точки использовался алгоритм “Лесенка Монтгомери”.<br/>

В GoogleClass(задание: отчет о вычислении кратных точек эллиптических кривых) лежит полный отчет о проделанной работе, теоритическая часть и скриншоты кода с пояснением.

## Сборка

Для запуска кода, у Вас должна быть установлена библиотека gmp, ссылка на скачивание библиотеки: https://gmplib.org.<br/>

## Тесты 

1 Проверить, что результирующая точка Q лежит на кривой.<br/>
2 Проверить, что [q]P = O, где q – порядок группы точек.<br/>
3 Проверить, что [q + 1]P = P и [q − 1]P = −P.<br/>
4 Для двух случайных k1, k2 проверить, что [k1]P + [k2]P = [k1 + k2]P.<br/>

![скриншот тестов](https://github.com/bulgvkov/hesse_curve/blob/main/screenshotOfTests.png)
