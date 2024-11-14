#ifndef __THREAD_SAFE_WRITER_HH_
#define __THREAD_SAFE_WRITER_HH_

#include <iostream>
#include <fstream>
#include <mutex>
#include <string>
#include <thread>

class ThreadSafeWriter {
public:
    ThreadSafeWriter(const std::string &filename);

    ~ThreadSafeWriter();

    void write(const std::string &message);

private:
    std::ofstream outFile;
    std::mutex fileMutex;
};

#endif /* __THREAD_SAFE_WRITER_HH_ */
