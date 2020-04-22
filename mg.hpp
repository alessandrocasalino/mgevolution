#include <iostream>
#include <string>

struct mg_cosmology
{
  std::vector<double> a;
  std::vector<double> H;
};

void mg_import(mg_cosmology * mg_cosmo) {

    ifstream fin;

    fin.open("mg_bk.csv", ios::in);

    int count = -1;

    vector<string> row;
    string line, word;

    do {

        row.clear();

        // read an entire row and
        // store it in a string variable 'line'
        getline(fin, line);

        // used for breaking words
        stringstream s(line);

        getline(s, word, ',');

        // read every column data of a row and
        // store it in a string variable, 'word'
        while (s >> word) {

            //cout << word << " ";
            row.push_back(word);
        }

        //cout << "\n";

        if(count>=0){

            mg_cosmo->a.push_back(stod(row[1]));
            mg_cosmo->H.push_back(stod(row[2]));

        }

        count++;

    } while (line != "");

}
