import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

machine_eps = np.finfo(float).eps

def get_forcing_Fourier_coefs(k, n, A=1, sigma=0.6): 
    return A * np.random.normal(size=n) * k**4 *np.exp(-0.5*sigma**2 *k**2)


def generate_forcing(samples, n, A=1):
    f = np.empty((samples, n))
    f_hat = np.empty((samples, 17))
    
    for k in range(0,17):
        f_hat[:, k] = get_forcing_Fourier_coefs(k, samples, A, sigma=0.6)

    # print(f_hat)
    
    # f_hat = np.array([get_forcing_Fourier_coefs(k, samples, A=0.001, sigma=0.6) for k in range(0,n)])

    f[:] = np.fft.ifft(f_hat, n=n).real

    # pbar.close()
    return f




# def generate_forcing(samples, n):
#     f = np.empty((samples, n))
#     i=1
#     failures = 0
#     # pbar = tqdm(total=samples)

#     while i < samples:
        
#         f_hat = [get_forcing_Fourier_coefs(k, A=0.001) for k in range(0,n)]
#         # print(f_hat[:17]) #only first 16 wavenumbers should be above machine precision
#         try:
#             assert (f_hat[16] < machine_eps and f_hat[15] > machine_eps)
#             i+=1
#             # pbar.update(1)
#         except AssertionError:
#             # failures += 1
#             # if failures > samples*100 and i==0:
#             #     raise RuntimeError("Could not find forcing in suitable time, modify sigma.")
#             #method to modify sigma if not finding suitable forcings?
#             continue

#         f[i-1] = np.fft.ifft(f_hat).real

#     # pbar.close()
#     return f





if __name__=="__main__":
    print("Machine precision is: ", machine_eps)
    n_points = 2048
    n_samples = 5

    f = generate_forcing(n_samples, n_points)

    for i in range(n_samples):
        plt.plot(np.linspace(0,1,n_points), f[i,:])

    plt.title("Spatial forcing from 16 Fourier modes - adapted from LaBryer et al.")
    plt.grid()
    plt.xlim(0,1)
    plt.show()