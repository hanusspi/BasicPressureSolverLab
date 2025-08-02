#pragma once
#include <iostream>

// Abstract Base Class (Interface) for all exercises
class Exercise {
public:
    virtual void run() = 0;  // Pure virtual function
    virtual ~Exercise() {}       // Virtual destructor for proper cleanup
};