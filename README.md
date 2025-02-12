
# Самописный РС кодер/декодер

## Использование

### Скрипт с примерами

```bash
make main
./main
```

### Скрипт с кодером

```bash
make coder
./coder
```

### Удаление исполняемых файлов

```bash
make clean
```

## Список задач

- [ ] **Рефакторинг существующих функций**:
  - [ ] Определиться с *gf_elem_t* и *gf_inner_t* (убрать или оставить)
  - [ ] Оптимизация функций по памяти (избавиться от лишних инициализаций)

- [ ] **Операции с элементами поля Галуа**:
  - [ ] Вычитание элементов (`sub`)
  - [ ] Деление элементов (`div`)
  - [ ] Возведение в целочисленную степень (`pow`)

- [ ] **Операции с полиномами**:
  - [ ] Подставление аргумента и вычисление значения (`call`)
  - [ ] Сложение полиномов (`add`)
  - [ ] Вычистение полиномов (`sub`)
  - [ ] Деления полиномов с остатком (`div`)

- [ ] **Прочее**:
  - [ ] Добавить тестирование существующих модулей
  - [ ] Оформить описание и документацию к библиотеке
  - [ ] Написать примеры использования
