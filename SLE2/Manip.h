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
    };        // �����, ������� ������������ ��� ������ � �������

    vector<string> columns;          // �������� ��������
    vector<vector<vector<string>>> rows;  // ��������� ������  [����� ������][����� �������][����� �����]
    vector<int> width;            // ������ ������� �������
    int lineLength = 0;            // ����� ������ (� ������ �������� '|' )
    void stable() {    // ������ �� ��������!!! ���������� ��� ����� ���������� ������
        for (int i = 0; i < rows.back()[0].size(); i++) {
            stringstream bufer(rows.back()[0][i]);  // ���� �������
            string line;
            int n = 0;
            while (getline(bufer, line)) {          // ��������� ������� ��������� �� ������ ������(�������� ��������� �������)
                n++;
                if (n > rows.back().size()) {    // ���� ������ ������ �� ������������� ���������� ���������� ����, �� ����������� ������ ������
                    rows.back().push_back(rows.back()[0]);  // ��������� ��� ���� "������" � ������
                    for (int j = 0; j < rows.back().back().size(); j++) {
                        rows.back().back()[j] = "";      // ���������� ����������� "������"
                    }
                }
                rows.back()[n - 1][i] = line;
            }
        }
    }
public:
    Writer(vector<string> columns)  // ������������� ��������
        : columns(columns)
        , width(vector<int>(columns.size()))
    {
        lineLength = 4 * columns.size() + 1;
        add(columns);
    }

    void add(vector<string> data) {  // ���������� ������ � �������
        if (data.size() != columns.size()) {
            throw exception("���������� ��������� �� ������������� ���������� ��������!");
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
                    string space(width[col] + 1 - (int)rows[line][lvl][col].size(), ' ');  // ���������� ��������, ������� ���������� ������� � �����
                    printf("|  %s%s", rows[line][lvl][col].c_str(), space.c_str());
                }
                printf("|\n");                              // ��������� ������� � ������� �� ����� ������
            }
            string lowerLine(lineLength, '_');                    // ����������� ������ (��������� ��� ������ ������)
            printf("\x1b[38;2;70;188;255m%s\n", lowerLine.c_str());

        }
        printf("\x1b[38;2;226;226;226m");
    }
};
