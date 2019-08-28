
#ifndef FrameLib_THREADING_H
#define FrameLib_THREADING_H

#include "FrameLib_Types.h"

#include <atomic>
#include <algorithm>
#include <chrono>
#include <thread>
#include <vector>

/**
 
 @defgroup Threading
 
 */

#ifdef __linux__

// Linux specific definitions

#include <pthread.h>
#include <semaphore.h>

namespace OS_Specific
{
    typedef pthread_t OSThreadType;
    typedef sem_t OSSemaphoreType;
    typedef void *OSThreadFunctionType(void *arg);
}

#elif defined(__APPLE__)

// OSX specific definitions

#include <pthread.h>
#include <mach/mach.h>

namespace OS_Specific
{
    typedef pthread_t OSThreadType;
    typedef semaphore_t OSSemaphoreType;
    typedef void *OSThreadFunctionType(void *arg);
}

#else
    
// Windows OS specific definitions
    
#include <windows.h>

namespace OS_Specific
{
    typedef HANDLE OSThreadType;
    typedef HANDLE OSSemaphoreType;
    typedef DWORD WINAPI OSThreadFunctionType(LPVOID arg);
}

#endif

// Helpers for exchange fixed values with atomic types

template <class T>
bool compareAndSwap(std::atomic<T>& value, T comparand, T exchange)
{
    return value.compare_exchange_strong(comparand, exchange);
}

template <class T>
bool nullSwap(std::atomic<T *>& value, T *exchange)
{
    T *comparand = nullptr;
    return value.compare_exchange_strong(comparand, exchange);
}


/**
 
 @class FrameLib_CountedPointer
 
 @ingroup Threading
 
 @brief a templated counted pointer for the purpose of lock-free operations.
 
 */

template <class T>
struct FrameLib_CountedPointer
{
    FrameLib_CountedPointer() : mPointer(nullptr), mCount(0) {}
    FrameLib_CountedPointer(T *item, uintptr_t count) : mPointer(item), mCount(count) {}
    
    bool operator==(const FrameLib_CountedPointer& a) const
    {
        return a.mPointer == mPointer && a.mCount == mCount;
    }
    
    // This structure must not have padding bits, so we use a pointer sized mCount
    
    T *mPointer;
    uintptr_t mCount;
};


/**
 
 @class FrameLib_LockFreePointer
 
 @ingroup Threading
 
 @brief a templated lock-free pointer that can be swapped safely.
 
 */

template <class T>
struct FrameLib_LockFreePointer : public std::atomic<FrameLib_CountedPointer<T>>
{
    using Pointer = FrameLib_CountedPointer<T>;
    using Base = std::atomic<Pointer>;
    
    FrameLib_LockFreePointer() : Base(Pointer()) {}
    
    FrameLib_LockFreePointer(const Pointer& pointer) : Base(pointer) {}
    
    bool trySwap(T *item, Pointer compare)
    {
        return compareAndSwap(*this, compare, Pointer(item, compare.mCount + 1));
    }
};


/**
 
 @class FrameLib_SpinLock
 
 @ingroup Threading

 @brief a spinlock that can be locked, attempted or acquired.
 
 */

class FrameLib_SpinLock
{
    
public:
    
    FrameLib_SpinLock() : mAtomicLock(false) {}
    ~FrameLib_SpinLock() { acquire(); }
    
    // Non-copyable
    
    FrameLib_SpinLock(const FrameLib_SpinLock&) = delete;
    FrameLib_SpinLock& operator=(const FrameLib_SpinLock&) = delete;
    
    bool attempt() { return compareAndSwap(mAtomicLock, false, true); }
    void acquire() { while (attempt() == false); }
    void release() { compareAndSwap(mAtomicLock, true, false); }
    
private:
    
    std::atomic<bool> mAtomicLock;
};


/**
 
 @class FrameLib_SpinLockHolder
 
 @ingroup Threading

 @brief a RAII hold utility for a FrameLib_SpinLock
 
 */

class FrameLib_SpinLockHolder
{
    
public:
    
    FrameLib_SpinLockHolder(FrameLib_SpinLock *lock) : mLock(lock) { if (mLock) mLock->acquire(); }
    ~FrameLib_SpinLockHolder() { if (mLock) mLock->release(); }
    
    // Non-copyable
    
    FrameLib_SpinLockHolder(const FrameLib_SpinLockHolder&) = delete;
    FrameLib_SpinLockHolder& operator=(const FrameLib_SpinLockHolder&) = delete;
    
    void destroy()
    {
        if (mLock)
            mLock->release();
        mLock = nullptr;
    }
    
private:
    
    FrameLib_SpinLock *mLock;
};


/**
 
 @class FrameLib_Thread
 
 @ingroup Threading

 @brief lightweight joinable thread with variable priority level
 
 The thread must be joined before destruction.
 
 */

class FrameLib_Thread
{
    typedef void ThreadFunctionType(void *);
    
public:
    
    enum PriorityLevel {kLowPriority, kMediumPriority, kHighPriority, kAudioPriority};

    FrameLib_Thread(PriorityLevel priority, ThreadFunctionType *threadFunction, void *arg)
    : mInternal(nullptr), mPriority(priority), mThreadFunction(threadFunction), mArg(arg), mValid(false)
    {}

    ~FrameLib_Thread();

    static unsigned int maxThreads()
    {
        return std::max(1U, std::thread::hardware_concurrency());
    }
    
    static void sleepCurrentThread(unsigned long nanoseconds)
    {
        std::this_thread::sleep_for(std::chrono::nanoseconds(nanoseconds));
    }
    
    static int currentThreadPriority();
    
    // Non-copyable
    
    FrameLib_Thread(const FrameLib_Thread&) = delete;
    FrameLib_Thread& operator=(const FrameLib_Thread&) = delete;
    
    void start();
    void join();

private:
    
    // threadStart is a quick OS-style wrapper to call the object which calls the relevant static function
    
    static OS_Specific::OSThreadFunctionType threadStart;
    void call() { mThreadFunction(mArg); }
    
    // Data
    
    OS_Specific::OSThreadType mInternal;
    PriorityLevel mPriority;
    ThreadFunctionType *mThreadFunction;
    void *mArg;
    bool mValid;
};


/**
 
 @class FrameLib_Semaphore
 
 @brief a semaphore class wrapping an OS-level semaphore
 
 @ingroup Threading

 The semaphore must be clsed before destruction.
 
 */

class FrameLib_Semaphore
{

public:
    
    FrameLib_Semaphore(long maxCount);
    ~FrameLib_Semaphore();
    
    // Non-copyable

    FrameLib_Semaphore(const FrameLib_Semaphore&) = delete;
    FrameLib_Semaphore& operator=(const FrameLib_Semaphore&) = delete;
    
    void close();
    void signal(long n);
    bool wait();
    
private:
    
    // Data
    
    OS_Specific::OSSemaphoreType mInternal;
    bool mValid;
};


/**
 
 @class FrameLib_TriggerableThread

 @ingroup Threading

 @brief a thread that can be triggered from another thread (there is no mechanism to check progress)
 
 The thread should be joined before destruction.
 
 */

class FrameLib_TriggerableThread
{
    
public:

    FrameLib_TriggerableThread(FrameLib_Thread::PriorityLevel priority) : mThread(priority, threadEntry, this), mSemaphore(1) {}
    virtual ~FrameLib_TriggerableThread() {}
    
    // Non-copyable
    
    FrameLib_TriggerableThread(const FrameLib_TriggerableThread&) = delete;
    FrameLib_TriggerableThread& operator=(const FrameLib_TriggerableThread&) = delete;
    
    // Start and join
    
    void start() { mThread.start(); }
    void join();
    
    // Trigger the thread to do something

    void signal() { mSemaphore.signal(1); };
    
private:
    
    // threadEntry simply calls threadClassEntry which calls the task handler
    
    static void threadEntry(void *thread);
    void threadClassEntry();
    
    // Override this and provide code for the thread's functionality
    
    virtual void doTask() = 0;
    
    // Data
    
    FrameLib_Thread mThread;
    FrameLib_Semaphore mSemaphore;
};


/**
 
 @class FrameLib_DelegateThread
 
 @ingroup Threading

 @brief a thread to delegate tasks to, which can be then be checked for completion
 
 The thread should be joined before desctruction.
 
 */

class FrameLib_DelegateThread
{
    
public:

    FrameLib_DelegateThread(FrameLib_Thread::PriorityLevel priority)
    : mThread(priority, threadEntry, this), mSemaphore(1), mSignaled(false), mFlag(0) {}
    virtual ~FrameLib_DelegateThread() {}

    // Non-copyable
    
    FrameLib_DelegateThread(const FrameLib_DelegateThread&) = delete;
    FrameLib_DelegateThread& operator=(const FrameLib_DelegateThread&) = delete;
    
    // Start and join

    void start() { mThread.start(); }
    void join();
    
    // Signal the thread to do something if it is not busy (returns true if the thread was signalled, false if busy)
    
    bool signal();
    
    // Check if the task has completed (false if the task has not been signalled or has not completed)
    
    bool completed();
    
     // Wait for thread's completion of a task - returns true  first time after a signal / false for subsequent calls / the thread has not been signalled
    
    bool waitForCompletion();
    
private:
    
    // threadEntry simply calls threadClassEntry which calls the task handler
    
    static void threadEntry(void *thread);
    void threadClassEntry();
    
    // Override this and provide code for the thread's functionality

    virtual void doTask() = 0;
    
    // Data

    FrameLib_Thread mThread;
    FrameLib_Semaphore mSemaphore;
    
    bool mSignaled;
    std::atomic<int> mFlag;
};


/**
 
 @class FrameLib_TriggerableThreadSet
 
 @ingroup Threading
 
 @brief a set of indexed thread to use as a worker pool
 
 The threads should be joined before desctruction.
 
 */

class FrameLib_TriggerableThreadSet
{
    struct IndexedThread : public FrameLib_Thread
    {
        IndexedThread(PriorityLevel priority, FrameLib_TriggerableThreadSet *owner, unsigned int index)
        : FrameLib_Thread(priority, threadEntry, this), mOwner(owner), mIndex(index) {}
        
        FrameLib_TriggerableThreadSet *mOwner;
        unsigned int mIndex;
    };
    
public:
    
    FrameLib_TriggerableThreadSet(FrameLib_Thread::PriorityLevel priority, unsigned int size);
    virtual ~FrameLib_TriggerableThreadSet() {}
    
    // Start and join
    
    void start();
    void join();
    
    // Trigger the threads to do something
    
    void signal(int32_t n) { mSemaphore.signal(n); };

    // Get the size
    
    unsigned int size() { return static_cast<unsigned int>(mThreads.size()); };

private:
    
    // Deleted
    
    FrameLib_TriggerableThreadSet(const FrameLib_TriggerableThreadSet&);
    FrameLib_TriggerableThreadSet& operator=(const FrameLib_TriggerableThreadSet&);
    
    // threadEntry simply calls threadClassEntry which calls the task handler
    
    static void threadEntry(void *thread);
    void threadClassEntry(unsigned int index);
    
    // Override this and provide code for the thread's functionality
    
    virtual void doTask(unsigned int index) = 0;
    
    // Data
    
    FrameLib_OwnedList<IndexedThread> mThreads;
    FrameLib_Semaphore mSemaphore;
};

#endif
