import os
def readFa(fa):
    '''
    @msg: 读取一个fasta文件
    @param fa {str}  fasta 文件路径
    @return: {generator} 返回一个生成器，能迭代得到fasta文件的每一个序列名和序列
    '''
    with open(fa,'r') as FA:
        seqName,seq='',''
        while 1:
            line=FA.readline()
            line=line.strip('\n')
            if (line.startswith('>') or not line) and seqName:
                yield((seqName,seq))
            if line.startswith('>'):
                seqName = line[1:]
                seq=''
            else:
                seq+=line
            if not line:break
    

def getSeq(fa,querySeqName,start=1,end=0):
    '''
    @msg: 获取fasta文件的某一条序列
    @param fa {str}  fasta 文件路径
    @param querySeqName {str}  序列名
    @param start {int}  截取该序列时，起始位置，可省略，默认为1
    @param end {int}  fasta 截取该序列时，最后位置，可省略，默认为该序列全长
    @return: {str} 返回找到(截取到)的序列
    '''
    if start<0: start=start+1
    for seqName,seq in readFa(fa):
        if querySeqName==seqName:
            if end!=0: returnSeq = seq[start-1:end];print(start-1)
            else: returnSeq = seq[start-1:]
            return returnSeq


def getReverseComplement(sequence):
    '''
    @msg: 获取反向互补序列
    @param sequence {str}  一段DNA序列
    @return: {str} 返回反向互补序列
    '''
    sequence = sequence.upper()
    sequence = sequence.replace('A', 't')
    sequence = sequence.replace('T', 'a')
    sequence = sequence.replace('C', 'g')
    sequence = sequence.replace('G', 'c')
    return sequence.upper()[::-1]

def getGC(sequence):
    '''
    @msg: 获取某一条序列的GC含量
    @param sequence {str}  一段DNA序列
    @return: {float} 返回GC含量
    '''
    sequence=sequence.upper()
    return (sequence.count("G")+sequence.count("C"))/len(sequence)


def readSeqByWindow(sequence,winSize,stepSize):
    '''
    @msg: 滑窗读取某一条序列
    @param sequence {str}  一段DNA序列
    @param winSize {int}  窗口大小
    @param stepSize {int}  步长
    @return: {generator}  返回一个生成器，可迭代得到该序列的每一个窗口序列
    '''
    if stepSize<=0: return False
    now = 0
    seqLen = len(sequence)
    while(now+winSize-stepSize<seqLen):
        yield sequence[now:now+winSize]
        now+=stepSize

def getGapPos(sequence):
    '''
    @msg: 获取某条序列中gap的位置
    @param sequence {str}  一段DNA序列
    @return: {list}  返回一个列表，列表中每个元素为每个gap的起始和结束位置
    '''
    Ns = {'N', 'n'}
    result = []
    i = 0
    for base in sequence:
        i += 1
        if not base in Ns: continue
        if len(result) == 0 : result.append([i,i])
        elif i - result[-1][1] == 1: result[-1][1] = i
        else: result.append([i,i])
    return result

if __name__ == '__main__':
    fa_path = 'C:\\Users\\Jeff\\OneDrive\\桌面\\genome.fa'
    os.chdir('C:\\Users\\Jeff\\OneDrive\\桌面\\')
    for seqName,seq in readFa(fa_path):
        seqLen = len(seq)
        GC = getGC(seq)
        with open("genome_addN.fa",'a+`') as newFA:
            newFA.write('>'+seqName+"\n")
            newFA.write('N'*2000+seq+'N'*2000+"\n")