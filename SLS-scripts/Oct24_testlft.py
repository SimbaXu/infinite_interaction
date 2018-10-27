import control as co
import scipy.io as io

data = io.loadmat('Oct24_testlft_data.mat')

P = co.ss(data['A'], data['B'], data['C'], data['D'])
K = co.ss(data['Ak'], data['Bk'], data['Ck'], data['Dk'])
pk = co.ss(data['Alft'], data['Blft'], data['Clft'], data['Dlft'])

pk_co = P.lft(K)
mag, phase, _ = pk_co.freqresp([1])
import IPython
if IPython.get_ipython() is None:
    IPython.embed()
