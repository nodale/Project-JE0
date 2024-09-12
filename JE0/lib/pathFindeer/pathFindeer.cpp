#include "pathFindeer.h"

std::filesystem::path getCurrentPath()
{
    std::filesystem::path currentPath = std::filesystem::current_path();
    return currentPath;
}