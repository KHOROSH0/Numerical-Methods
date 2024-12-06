#pragma once
#pragma once

#include <vector>
#include <string>
#include <clocale>
#include <sstream>
#include <iostream>
using namespace std;

class Writer
{
private:
    vector<string> color{
      "\x1b[38;2;226;226;226m",
      "\x1b[38;2;70;188;255m",
    };        // цвета, которые используются при выводе в консоль

    vector<string> columns;          // название столбцов
    vector<vector<vector<string>>> rows;  // выводимые строки  [номер строки][номер столбца][номер этажа]
    vector<int> width;            // ширина каждого столбца
    int lineLength = 0;            // длина строки (с учётом символов '|' )
    void stable() {    // самому не вызывать!!! Вызывается сам после добавления строки
        for (int i = 0; i < rows.back()[0].size(); i++) {
            stringstream bufer(rows.back()[0][i]);  // берём столбец
            string line;
            int n = 0;
            while (getline(bufer, line)) {          // проверяем сколько переводов на другую строку(выделяем отдельный столбец)
                n++;
                if (n > rows.back().size()) {    // если размер строки не соответствует количество разделённых слов, то увеличиваем размер строки
                    rows.back().push_back(rows.back()[0]);  // добавляем ещё одну "ячейку" в строке
                    for (int j = 0; j < rows.back().back().size(); j++) {
                        rows.back().back()[j] = "";      // опустошаем добавленную "ячейку"
                    }
                }
                rows.back()[n - 1][i] = line;
            }
        }
    }
public:
    Writer(vector<string> columns)  // инициализация столбцов
        : columns(columns)
        , width(vector<int>(columns.size()))
    {
        lineLength = 4 * columns.size() + 1;
        add(columns);
    }

    void add(vector<string> data) {  // добавление строки в таблицу
        if (data.size() != columns.size()) {
            throw exception("Количество элементов не соответствует количеству столбцов!");
        }
        for (int i = 0; i < columns.size(); i++) {
            stringstream line(data[i]);
            string word;
            while (getline(line, word)) {
                if ((int)word.size() > width[i]) {
                    lineLength += (int)word.size() - width[i];
                    width[i] = (int)word.size();
                }
            }
        }
        rows.push_back({ data });
        stable();
    }

    void print() {
        setlocale(LC_ALL, "ru");
        for (int line = 0; line < rows.size(); line++) {
            for (int lvl = 0; lvl < rows[line].size(); lvl++) {
                printf(color[line % color.size()].c_str());
                for (int col = 0; col < rows[line][lvl].size(); col++) {
                    string space(width[col] + 1 - (int)rows[line][lvl][col].size(), ' ');  // количество пробелов, которые необходимо вывести в конце
                    printf("|  %s%s", rows[line][lvl][col].c_str(), space.c_str());
                }
                printf("|\n");                              // закрываем столбец и переход на новую строку
            }
            string lowerLine(lineLength, '_');                    // разделяющая строка (разделяет две разыне строки)
            printf("\x1b[38;2;70;188;255m%s\n", lowerLine.c_str());

        }
        printf("\x1b[38;2;226;226;226m");
    }
};
