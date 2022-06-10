#include "dimacs.h"


DIMACS::DIMACS(const std::string filename)
{
    std::ifstream file;
    file.open(filename);

    if(file.fail()){
        exit(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(file, line))
    {
        if (line[0] == 'c')
        {
        }
        else if (line[0] == 'p')
        {
            std::istringstream iss(line);
            std::string element;
            size_t counter{0};
            while (std::getline(iss, element, ' '))
            {
                if (counter == 2)
                {
                    std::istringstream(element) >> vertexCount;
                }
                else if (counter == 3)
                {
                    std::istringstream(element) >> edgeCount;
                }

                ++counter;
            }

            matrix = std::make_unique<Matrix<bool>>(vertexCount, vertexCount, false);
        }
        else if (line[0] == 'e')
        {
            size_t x;
            size_t y;

            std::istringstream iss(line);
            std::string element;
            size_t counter{0};
            while (std::getline(iss, element, ' '))
            {
                if (counter == 1)
                {
                    std::istringstream(element) >> x;
                }
                else if (counter == 2)
                {
                    std::istringstream(element) >> y;
                }

                ++counter;
            }

            (*matrix)(x - 1, y - 1) = true;
            (*matrix)(y - 1, x - 1) = true;
        }
    }

    file.close();
}
