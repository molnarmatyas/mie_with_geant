#include "ThreadSafeWriter.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
namespace { G4Mutex myMutex = G4MUTEX_INITIALIZER; }

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
  G4AutoLock lock(&myMutex);
	if (outFile.is_open())
	{
		outFile.close();
	}
}

void ThreadSafeWriter::write(const std::string &message)
{
	//std::lock_guard<std::mutex> guard(fileMutex);
  G4AutoLock lock(&myMutex);
	outFile << message << std::endl;
}
