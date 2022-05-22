import numpy as np
import matplotlib.pyplot as plt
import pandas

# Figs and plots settings
SMALL_SIZE = 14
BIG_SIZE = 18

plt.rc('font', size=BIG_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIG_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=BIG_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=BIG_SIZE)    # legend fontsize
plt.rc('figure', titlesize=BIG_SIZE)  # fontsize of the figure title

# # Strain values of simulations
# kefile = open("data/results/kin_eng_dt01.txt", 'r')
# time_dt01 = []
# ke_dt01 = []
# for line in kefile.readlines():
#     time_dt01.append(float(line.split(',')[0]))
#     ke_dt01.append(float(line.split(',')[1]))

# kefile.close() 

# kefile = open("data/results/kin_eng_dt001.txt", 'r')
# time_dt001 = []
# ke_dt001 = []
# for line in kefile.readlines():
#     time_dt001.append(float(line.split(',')[0]))
#     ke_dt001.append(float(line.split(',')[1]))

# kefile.close() 

# kefile = open("data/results/kin_eng_dt0001.txt", 'r')
# time_dt0001 = []
# ke_dt0001 = []
# for line in kefile.readlines():
#     time_dt0001.append(float(line.split(',')[0]))
#     ke_dt0001.append(float(line.split(',')[1]))

# kefile.close() 

# kefile = open("data/results/kin_eng_sincor_dt001.txt", 'r')
# time_sc_dt001 = []
# ke_sc_dt001 = []
# for line in kefile.readlines():
#     time_sc_dt001.append(float(line.split(',')[0]))
#     ke_sc_dt001.append(float(line.split(',')[1]))

# kefile.close() 

# plt.figure(figsize=(11, 4), dpi=200)
# plt.plot(time_dt01, ke_dt01, 'r', time_dt001, ke_dt001, 'b', time_dt0001, ke_dt0001, 'k')
# plt.tight_layout()
# plt.title("Kinetic Energy of the system")
# plt.legend(["dt = 0.01 sec", "dt = 0.001 sec", "dt = 0.0001 sec"])
# plt.xlabel("Time (s)")
# plt.ylabel("Kinetic Energy (J)")
# # plt.ylim(-1e-7, 7e-7)
# # plt.xlim(2251.5, 2263)
# # plt.show() 
# plt.savefig('kin_eng_dt_study.png')
# plt.close()

# plt.figure(figsize=(11, 4), dpi=200)
# plt.plot(time_dt0001, ke_dt0001, 'k')
# plt.tight_layout()
# plt.title("Kinetic Energy of the system")
# plt.legend(["dt = 0.0001 sec"])
# plt.xlabel("Time (s)")
# plt.ylabel("Kinetic Energy (J)")
# # plt.ylim(-1e-7, 7e-7)
# # plt.xlim(2251.5, 2263)
# # plt.show() 
# plt.savefig('kin_eng_dt0001.png')
# plt.close()

# plt.figure(figsize=(11, 4), dpi=200)
# plt.plot(time_sc_dt001, ke_sc_dt001, 'k')
# plt.tight_layout()
# plt.title("Kinetic Energy of the system")
# plt.legend(["dt = 0.001 sec Single-Core"])
# plt.xlabel("Time (s)")
# plt.ylabel("Kinetic Energy (J)")
# # plt.ylim(-1e-7, 7e-7)
# # plt.xlim(2251.5, 2263)
# # plt.show() 
# plt.savefig('kin_eng__sc_dt001.png')
# plt.close()

# kefile = open("data/results/kin_eng.txt", 'r')
# time_cc10t10dt001 = []
# ke_cc10t10dt001 = []
# for line in kefile.readlines():
#     time_cc10t10dt001.append(float(line.split(',')[0]))
#     ke_cc10t10dt001.append(float(line.split(',')[1]))

# kefile.close() 

# plt.figure(figsize=(11, 4), dpi=200)
# plt.plot(time_cc10t10dt001, ke_cc10t10dt001, 'k')
# plt.tight_layout()
# plt.title("Kinetic Energy: cc10 t10 dt001 f2")
# # plt.legend(["dt = 0.001 sec Single-Core"])
# plt.xlabel("Time (s)")
# plt.ylabel("Kinetic Energy (J)")
# # plt.ylim(-1e-7, 7e-7)
# # plt.xlim(2251.5, 2263)
# # plt.show() 
# plt.savefig('kin_eng_cc10t10dt001f2.png')
# plt.close()

# kefile = open("data/kin_eng_tests/kin_engcc10t10dt01f2.txt", 'r')
# time_cc10t10dt01f2 = []
# ke_cc10t10dt01f2 = []
# for line in kefile.readlines():
#     time_cc10t10dt01f2.append(float(line.split(',')[0]))
#     ke_cc10t10dt01f2.append(float(line.split(',')[1]))

# kefile.close() 

# kefile = open("data/kin_eng_tests/kin_engcc10t10dt005f2.txt", 'r')
# time_cc10t10dt005f2 = []
# ke_cc10t10dt005f2 = []
# for line in kefile.readlines():
#     time_cc10t10dt005f2.append(float(line.split(',')[0]))
#     ke_cc10t10dt005f2.append(float(line.split(',')[1]))

# kefile.close() 

# kefile = open("data/kin_eng_tests/kin_engcc10t10dt002f2.txt", 'r')
# time_cc10t10dt002f2 = []
# ke_cc10t10dt002f2 = []
# for line in kefile.readlines():
#     time_cc10t10dt002f2.append(float(line.split(',')[0]))
#     ke_cc10t10dt002f2.append(float(line.split(',')[1]))

# kefile.close() 

# kefile = open("data/kin_eng_tests/kin_engcc10t10dt001f2.txt", 'r')
# time_cc10t10dt001f2 = []
# ke_cc10t10dt001f2 = []
# for line in kefile.readlines():
#     time_cc10t10dt001f2.append(float(line.split(',')[0]))
#     ke_cc10t10dt001f2.append(float(line.split(',')[1]))

# kefile.close() 

# kefile = open("data/kin_eng_tests/kin_engcc10t10dt0005f2.txt", 'r')
# time_cc10t10dt0005f2 = []
# ke_cc10t10dt0005f2 = []
# for line in kefile.readlines():
#     time_cc10t10dt0005f2.append(float(line.split(',')[0]))
#     ke_cc10t10dt0005f2.append(float(line.split(',')[1]))

# kefile.close() 

# kefile = open("data/kin_eng_tests/kin_engcc10t10dt0001f2.txt", 'r')
# time_cc10t10dt0001f2 = []
# ke_cc10t10dt0001f2 = []
# for line in kefile.readlines():
#     time_cc10t10dt0001f2.append(float(line.split(',')[0]))
#     ke_cc10t10dt0001f2.append(float(line.split(',')[1]))

# kefile.close() 

# plt.figure(figsize=(11, 4), dpi=200)
# plt.plot(time_cc10t10dt005f2, ke_cc10t10dt005f2, 'b',
#         time_cc10t10dt002f2, ke_cc10t10dt002f2, 'g', time_cc10t10dt001f2, ke_cc10t10dt001f2, 'r',
#         time_cc10t10dt0005f2, ke_cc10t10dt0005f2, 'y', time_cc10t10dt0001f2, ke_cc10t10dt0001f2, 'm')
# plt.tight_layout()
# plt.title("Kinetic Energy: cc10 t10 dt001 f2")
# plt.legend(["dt = 0.005 sec", "dt = 0.002 sec", "dt = 0.001 sec", "dt = 0.0005 sec", "dt = 0.0001 sec"])
# plt.xlabel("Time (s)")
# plt.ylabel("Kinetic Energy (J)")
# plt.ylim(0, 0.0001)
# plt.xlim(0, 2)
# # plt.show() 
# plt.savefig('kin_eng_cc10t10f2_dt_study.png')
# plt.close()





# kefile = open("data/results/kin_eng.txt", 'r')
# time = []
# ke = []
# for line in kefile.readlines():
#     time.append(float(line.split(',')[0]))
#     ke.append(float(line.split(',')[1]))

# kefile.close() 

# plt.figure(figsize=(11, 4), dpi=200)
# plt.plot(time, ke, 'b')
# plt.tight_layout()
# plt.title("Kinetic Energy of the system")
# # plt.legend(["dt = 0.01 sec", "dt = 0.001 sec", "dt = 0.0001 sec"])
# plt.xlabel("Time (s)")
# plt.ylabel("Kinetic Energy (J)")
# # plt.ylim(-1e-7, 7e-7)
# # plt.xlim(2251.5, 2263)
# # plt.show() 
# plt.savefig('kin_eng_250clips.png')
# plt.close()


# 10,   0.00001,  10,      15,     10361
# 10,   0.00001,  20,      15,     11026
# 10,   0.00001,  30,      15,     11842
# 10,   0.000005, 50,      15,     27728
# 10,   0.000002, 70,      15,     71474

nclips = [10, 20, 30, 50, 70]
simtime = [10361, 11026, 11842, 27728, 71474-0*3600]

plt.figure(figsize=(11, 4), dpi=200)
plt.plot(nclips, simtime, 'bs-')
plt.tight_layout()
plt.title("Hanging chain: f = 15Hz, T = 10 sec")
# plt.legend(["dt = 0.01 sec", "dt = 0.001 sec", "dt = 0.0001 sec"])
plt.xlabel("num of clips")
plt.ylabel("Run time (s)")
plt.text(8.0, 11500, 'dt = 1E-4', fontsize=12)
plt.text(18.0, 12500, 'dt = 1E-4', fontsize=12)
plt.text(28.0, 14000, 'dt = 1E-4', fontsize=12)
plt.text(45.0, 30000, 'dt = 5E-5', fontsize=12)
plt.text(63.0, 70000, 'dt = 2E-5', fontsize=12)
# plt.ylim(-1e-7, 7e-7)
# plt.xlim(2251.5, 2263)
# plt.show()
plt.savefig('multicore_runtime.png')
plt.close()