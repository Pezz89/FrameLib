
#ifndef FRAMELIB_FROMHOST_H
#define FRAMELIB_FROMHOST_H

#include "FrameLib_DSP.h"
#include "FrameLib_HostProxy.h"

class FrameLib_FromHost final : public FrameLib_Processor
{
    // Typedef for concision
    
    typedef std::unique_ptr<std::vector<double>> OwnedFrame;
    
    // A FIFO list for storing parameter frame additions
    
    struct SerialList
    {
        SerialList() : mTop(nullptr), mTail(nullptr) {}
        ~SerialList() { clear(); }
        
        struct Item
        {
            Item() : mNext(nullptr) {}
            Item(const FrameLib_Parameters::Serial& serial) : mSerial(serial), mNext(nullptr) {}
            Item(const SerialList& list) : mNext(nullptr)
            {
                for (Item *item = list.mTop; item; item = item->mNext)
                    mSerial.write(&item->mSerial);
            }
            
            FrameLib_Parameters::AutoSerial mSerial;
            Item *mNext;
        };
        
        void push(Item *item)
        {
            assert (item->mNext == nullptr && "item already in a list");
            
            if (mTail)
            {
                mTail->mNext = item;
                mTail = item;
            }
            else
                mTop = mTail = item;
        }
        
        Item *pop()
        {
            Item *item = mTop;
            
            mTop = mTop ? mTop->mNext : nullptr;
            mTail = (mTail == item) ? nullptr : mTail;
            
            if (item)
                item->mNext = nullptr;
            
            return item;
        }
        
        bool empty() const
        {
            return mTop == nullptr;
        }
        
        unsigned long size() const
        {
            unsigned long summedSize = 0;
            
            for (Item *item = mTop; item; item = item->mNext)
                summedSize += item->mSerial.size();
                
            return summedSize;
        }
        
        void clear()
        {
            for (Item *item = pop(); item; item = pop())
                delete item;
        }
        
        void reassign(SerialList& list)
        {
            if (!mTop)
                mTop = list.mTop;
            else
                mTail->mNext = list.mTop;
            if (list.mTail)
                mTail = list.mTail;
            
            list.mTop = nullptr;
            list.mTail = nullptr;
        }
        
    private:
        
        // Deleted

        SerialList(const SerialList&) = delete;
        SerialList& operator=(const SerialList&) = delete;
        
        // Data
        
        Item *mTop;
        Item *mTail;
    };
    
public:

    // The owner should inherit from this class and use these calls to send to all registered objects

    struct Proxy : public FrameLib_HostProxy<FrameLib_FromHost>
    {
        Proxy(bool copyStreams) : mCopyStreams(copyStreams) {}
        
        // Send a vector frame
        
        void sendFromHost(unsigned long index, const double *values, unsigned long N);
        void sendFromHost(unsigned long index, unsigned long stream, const double *values, unsigned long N);
        
        // Send a parameter frame
        
        void sendFromHost(unsigned long index, const FrameLib_Parameters::Serial *serial);
        void sendFromHost(unsigned long index, unsigned long stream, const FrameLib_Parameters::Serial *serial);
        
        // Send a parameter that takes a string
        
        void sendFromHost(unsigned long index, const char *tag, const char *string);
        void sendFromHost(unsigned long index, unsigned long stream, const char *tag, const char *string);
        
        // Send a parameter that takes a vector
        
        void sendFromHost(unsigned long index, const char *tag, const double *values, unsigned long N);
        void sendFromHost(unsigned long index, unsigned long stream, const char *tag, const double *values, unsigned long N);
        
        // Copy data from the first stream to another stream
        
        void copyData(void *streamOwner, unsigned long stream);
        
        bool mCopyStreams;
    };
    
private:
    
    // Parameter Info and Enums
    
    enum ParameterList { kMode };
    enum Modes { kValues, kParams };
    
    struct ParameterInfo : public FrameLib_Parameters::Info { ParameterInfo(); };

public:
    
    // Constructor / Destructor
    
    FrameLib_FromHost(FrameLib_Context context, FrameLib_Parameters::Serial *serialisedParameters, FrameLib_Proxy *proxy);
    ~FrameLib_FromHost();
    
    // Stream Awareness
    
    void setStream(void *streamOwner, unsigned long stream) override;
    
    // Info
    
    std::string objectInfo(bool verbose) override;
    std::string inputInfo(unsigned long idx, bool verbose) override;
    std::string outputInfo(unsigned long idx, bool verbose) override;

private:
    
    void process() override;
    
    // Swapping data with the proxy
    
    OwnedFrame swapVectorFrame(OwnedFrame& swapVector);
    void updateSerialFrame(SerialList &freeList, SerialList::Item *addSerial);

// Data
    
    FrameLib_SpinLock mLock;
    OwnedFrame mVectorFrame;
    SerialList mSerialFrame;
    SerialList mSerialFreeFrame;
    Modes mMode;
    
    Proxy *mProxy;
    void *mStreamOwner;
    unsigned long mStream;
    
    static ParameterInfo sParamInfo;
};

#endif
