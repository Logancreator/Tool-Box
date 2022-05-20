import io
import numpy as np
import logging

logging.basicConfig(filemode='log.txt', filename='fail.log', format="%(message)s")
logger = logging.getLogger(__name__)
'''
    getLogger()定义一个日志器
'''


class RefMatrix:
    """ 计算比对的概率 """

    def __init__(self, ref_seq_path):
        self.ref_seq_path = ref_seq_path
        ref_seq = self.read_seq(ref_seq_path)
        self.ref_length = len(ref_seq)
        self.init_matrix(ref_seq=ref_seq)

    def init_matrix(self, ref_seq):
        """ 根据参考序列的长度初始化两个向量
            - 第一个向量放ref
            - 第二个向量记录变异次数
        """
        self.ref, self.count = np.array(ref_seq, dtype=np.int8), np.zeros(
            self.ref_length, dtype=np.uint16)

    def read_seq(self, seq_path):
        """ 读取序列并转换成数字 """  ##pysam  bam格式录入
        with open(seq_path, 'rb') as fp:
            return bytearray().join(io.BytesIO(fp.read()))

    def safe_range(self, range):
        """ 计算之前，看一下范围是不是超过了ref的长度 """
        return range.start >= 0 and range.stop <= self.ref_length

    def check(self, reads, range, reads_id):
        """ 查找变异碱基，并在变异碱基的位置count值+1 """
        if self.safe_range(range):
            self.count[range] += ~self.ref[range] & reads == reads
        logger.error(reads_id)


    def save_result(self, result_path):
        """ 把计算完的结果保存下来 """
        import pickle
        with open(result_path, 'wb') as fp:
            fp.write(pickle.dumps(self.count))



if __name__ == '__main__':
    ref_seq_path = './sum/test.fa'
    ref_matrix = RefMatrix(ref_seq_path=ref_seq_path)
    ref_matrix.check(reads=[65, 10, 65], range=slice(0, 3), reads_id='xxxxxxx')
    ref_matrix.save_result('result.txt')


