#include "ThreadSafeWriter.hh"

ThreadSafeWriter::ThreadSafeWriter(const std::string &filename) : outFile(filename, std::ios::app)
{
	outFile = std::ofstream(filename, std::ios::app);
	if (!outFile.is_open())
	{
		throw std::ios_base::failure("Failed to open file");
	}
}

ThreadSafeWriter::~ThreadSafeWriter()
{
	if (outFile.is_open())
	{
		outFile.close();
	}
}

void ThreadSafeWriter::write(const std::string &message)
{
	std::lock_guard<std::mutex> guard(fileMutex);
	outFile << message << std::endl;
}