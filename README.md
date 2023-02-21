# gprMax_project

В проекте реализованы методы для решения обрадной задачи обработки георадарных данных с помощью нейронных сетей.
Презентация и результаты - Guide/ML_gpr.pptx

_________________________________________________________
ML.ipynb - обучение нейронки, описание внутри.
Используемые технологии: сети U-NET, AE, PSP-Net, simple_U-NET (взято из публикации), Seq_NET, CAE, VCAE, librosa-featues (фичи с либрозы)

Реализованно:
1. Решение обратной задачи
1.1 Тесты с автокодировщиками
2. Решение задачи повышения разрешение записи по горизонтали
3. Задача деконволюции (НЕ РАБОТАЕТ)

Известные проблемы и планы на будущее:
1. CAE - подать во внутреннее пространство информацию об проводимости (Т.е предположем, что мы знаем где находятся металлические объекты и хотим это учесть в предсказании)
2. VAE - приспособить эту архитектуру для поиска ТОЛЬКО проводников
3. Реализовать генетику на архитекруре simple_U-NET
4. Проблема с датасетом_2, может стоит отдельно его преобразовать и залить
5. Реальзовать lstm сети, аначале на простом датасете с изменением размерности временной записи
________________________________________________________
Датесеты находятся в облаке!!!
Dataset/ - сгенерированный датасет_1 - простая 2D модель
Dataset_2/ - сгенерированный датасет_2 - сложная 2D модель
Dataset_3/ - сгенерированный датасет_3
________________________________________________________
gprmax_library.py - библиотека с функциями и классов для моделирования (решения прямой задачи).
Известные проблемы:
1. при конвертации в numpy массив размерность вместо 128/128 реализует 126/128
2. пересмотреть смысл использования numpy массива, .png формат лучше для восприятия данных и публикации датасета
________________________________________________________
Tests/test_gprmax_library.py - тесты библиокеки по моделированию:
Известные проблемы:
1. Мало тестов
2. В тест с numpy-массивом добавить оценку размерности
________________________________________________________
generation_models.py - код генерации модели
Известные проблемы:
1. отсутствуют тесты
2. !!!!!!! Проблема с ошибкой к доступу памяти, решается конструкцией try-exception. Возможно на ubunty таких проблем не будет
________________________________________________________
gprMax_colab_template.ipynb - код для запуска gprMax в Colabе
________________________________________________________
terminal.txt - шпора для команд в терминале



