import random 
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from decimal import Decimal

Sex = np.zeros((16,16))
Yex = np.zeros((16,16))
SYc = np.zeros((16,16))

Y = Decimal('-0.6')
for i in range(16):
	Y += Decimal('0.1')
	S = Decimal('-0.6')
	for j in range(16):
		S += Decimal('0.1')
		Yext = 0
		Sext = 0
		YSco = 0
		Sfreq = []
		for x in range(10):			
			n = []
			Xm_f = []
			Ym_f = []
			Sm_f = []
			Xf_f = []
			Yf_f = []
			Sf_f = []
			S_f = []
			N_female = []
			N_male = []
			N_1YSinF = []
			N_2YSinF = []
			N_1YSinM = []
			N_2YSinM = []
			N_3YSinM = []
			
			population = 100
			generations = 100
			
			#Initial number of female
			Num_XX   = 0
			Num_XXY  = 18
			Num_XXS  = 15
			Num_XXYY = 3
			Num_XXSS = 6
			Num_XXYS = 6
			
			#Initial number of male
			Num_XY  = 12
			Num_XS  = 9
			Num_XYY = 8
			Num_XSS = 7
			Num_XYS = 13
			Num_XYYY = 0
			Num_XYYS = 0
			Num_XYSS = 1
			Num_XSSS = 0
			
			#Female fitness parameters
			s = float(S)
			y = float(Y)
			h = 0.50
			
			#Relative female fitness
			RF_XX   = 0.00
			RF_XXS  = 1.00-(h*s)
			RF_XXY  = 1.00-(h*y)
			RF_XXSS = 1.00-s
			RF_XXYS = 1.00-s-(h*(y-s))
			RF_XXYY = 1.00-y
			
			#Relative male fitness
			RF_XY   = 1.00
			RF_XS   = 0.00
			RF_XYY  = 1.00
			RF_XSS  = 0.00
			RF_XYS  = 1.00
			RF_XYYY = 0.00
			RF_XYYS = 0.10
			RF_XYSS = 0.10
			RF_XSSS = 0.00
			
			Total_number_of_chromosomes_in_female = (2*Num_XX)+(3*Num_XXY)+(3*Num_XXS)+(4*Num_XXYY)+(4*Num_XXSS)+(4*Num_XXYS)
			Total_number_of_chromosomes_in_male = (2*Num_XY)+(2*Num_XS)+(3*Num_XYY)+(3*Num_XSS)+(3*Num_XYS)+(4*Num_XYYY)+(4*Num_XYYS)+(4*Num_XYSS)+(4*Num_XSSS)
			Total_number_of_X_in_female = (2*Num_XX)+(2*Num_XXY)+(2*Num_XXS)+(2*Num_XXYY)+(2*Num_XXSS)+(2*Num_XXYS)
			Total_number_of_Y_in_female = Num_XXY+(2*Num_XXYY)+Num_XXYS
			Total_number_of_S_in_female = Num_XXS+(2*Num_XXSS)+Num_XXYS
			Total_number_of_X_in_male = Num_XY+Num_XS+Num_XYY+Num_XSS+Num_XYS+Num_XYYY+Num_XYYS+Num_XYSS+Num_XSSS
			Total_number_of_Y_in_male = Num_XY+(2*Num_XYY)+Num_XYS+(3*Num_XYYY)+(2*Num_XYYS)+Num_XYSS
			Total_number_of_S_in_male = Num_XS+(2*Num_XSS)+Num_XYS+Num_XYYS+(2*Num_XYSS)+(3*Num_XSSS)
			
			Xf_freq = Total_number_of_X_in_female/Total_number_of_chromosomes_in_female
			Yf_freq = Total_number_of_Y_in_female/Total_number_of_chromosomes_in_female
			Sf_freq = Total_number_of_S_in_female/Total_number_of_chromosomes_in_female
			Xm_freq = Total_number_of_X_in_male/Total_number_of_chromosomes_in_male
			Ym_freq = Total_number_of_Y_in_male/Total_number_of_chromosomes_in_male
			Sm_freq = Total_number_of_S_in_male/Total_number_of_chromosomes_in_male
			S_freq = (Total_number_of_S_in_female+Total_number_of_S_in_male)/(Total_number_of_S_in_female+Total_number_of_S_in_male+Total_number_of_Y_in_female+Total_number_of_Y_in_male)
			
			num_female = Num_XX+Num_XXY+Num_XXS+Num_XXYY+Num_XXSS+Num_XXYS
			num_male = Num_XY+Num_XS+Num_XYY+Num_XSS+Num_XYS+Num_XYYY+Num_XYYS+Num_XYSS+Num_XSSS
			
			egg_X = (RF_XX*6*Num_XX)+(RF_XXY*3*Num_XXY)+(RF_XXS*3*Num_XXS)
			egg_XY = (RF_XXY*3*Num_XXY)+(RF_XXYY*6*Num_XXYY)+(RF_XXYS*3*Num_XXYS)
			egg_XS = (RF_XXS*3*Num_XXS)+(RF_XXSS*6*Num_XXSS)+(RF_XXYS*3*Num_XXYS)
			Total_number_of_female_gametes = egg_X+egg_XY+egg_XS
			freq_egg_X = egg_X/Total_number_of_female_gametes
			freq_egg_XY = egg_XY/Total_number_of_female_gametes
			freq_egg_XS = egg_XS/Total_number_of_female_gametes
			
			sperm_X = (RF_XY*3*Num_XY)+(RF_XS*3*Num_XS)+(RF_XYY*1*Num_XYY)+(RF_XSS*1*Num_XSS)+(RF_XYS*1*Num_XYS)
			sperm_XY = (RF_XYY*2*Num_XYY)+(RF_XYS*1*Num_XYS)+(RF_XYYY*3*Num_XYYY)+(RF_XYYS*2*Num_XYYS)+(RF_XYSS*1*Num_XYSS)
			sperm_XS = (RF_XSS*2*Num_XSS)+(RF_XYS*1*Num_XYS)+(RF_XYYS*1*Num_XYYS)+(RF_XYSS*2*Num_XYSS)+(RF_XSSS*3*Num_XSSS)
			sperm_Y = (RF_XY*3*Num_XY)+(RF_XYY*2*Num_XYY)+(RF_XYS*1*Num_XYS)
			sperm_S = (RF_XS*3*Num_XS)+(RF_XSS*2*Num_XSS)+(RF_XYS*1*Num_XYS)
			sperm_YY = (RF_XYY*1*Num_XYY)+(RF_XYYY*3*Num_XYYY)+(RF_XYYS*1*Num_XYYS)
			sperm_SS = (RF_XSS*1*Num_XSS)+(RF_XYSS*1*Num_XYSS)+(RF_XSSS*3*Num_XSSS)
			sperm_YS = (RF_XYS*1*Num_XYS)+(RF_XYYS*2*Num_XYYS)+(RF_XYSS*2*Num_XYSS)
			Total_number_of_male_gametes = sperm_X+sperm_XY+sperm_XS+sperm_Y+sperm_S+sperm_YY+sperm_SS+sperm_YS
			freq_sperm_X = sperm_X/Total_number_of_male_gametes
			freq_sperm_XY = sperm_XY/Total_number_of_male_gametes
			freq_sperm_XS = sperm_XS/Total_number_of_male_gametes
			freq_sperm_Y = sperm_Y/Total_number_of_male_gametes
			freq_sperm_S = sperm_S/Total_number_of_male_gametes
			freq_sperm_YY = sperm_YY/Total_number_of_male_gametes
			freq_sperm_SS = sperm_SS/Total_number_of_male_gametes
			freq_sperm_YS = sperm_YS/Total_number_of_male_gametes
			
			Num_1YSinF = Num_XXY+Num_XXS
			Num_2YSinF = Num_XXYY+Num_XXSS+Num_XXYS
			Num_1YSinM = Num_XY+Num_XS
			Num_2YSinM = Num_XYY+Num_XSS+Num_XYS
			Num_3YSinM = Num_XYYY+Num_XYYS+Num_XYSS+Num_XSSS
			
			n.append(0)
			N_female.append(num_female)
			N_male.append(num_male)
			N_1YSinF.append(Num_1YSinF)
			N_2YSinF.append(Num_2YSinF)
			N_1YSinM.append(Num_1YSinM)
			N_2YSinM.append(Num_2YSinM)
			N_3YSinM.append(Num_3YSinM)
			Xf_f.append(Xf_freq)
			Yf_f.append(Yf_freq)
			Sf_f.append(Sf_freq)
			Xm_f.append(Xm_freq)
			Ym_f.append(Ym_freq)
			Sm_f.append(Sm_freq)
			S_f.append(S_freq)
			
			for x in range(generations):
				n.append(x+1)
				num_female = 0
				num_male = 0
				Num_XX   = 0
				Num_XXY  = 0
				Num_XXS  = 0
				Num_XXYY = 0
				Num_XXSS = 0
				Num_XXYS = 0
				Num_XY  = 0
				Num_XS  = 0
				Num_XYY = 0
				Num_XSS = 0
				Num_XYS = 0
				Num_XYYY = 0
				Num_XYYS = 0
				Num_XYSS = 0
				Num_XSSS = 0
				egg_X = 0
				egg_XY = 0
				egg_XS = 0
				sperm_X = 0
				sperm_XY = 0
				sperm_XS = 0
				sperm_Y = 0
				sperm_S = 0
				sperm_YY = 0
				sperm_SS = 0
				sperm_YS = 0
				
				for x in range(population):
					egg = np.random.choice([10, 20, 30], 1, p=[freq_egg_X, freq_egg_XY, freq_egg_XS])
					sperm = np.random.choice([1, 2, 3, 4, 5, 6, 7, 8], 1, p=[freq_sperm_X, freq_sperm_XY, freq_sperm_XS, freq_sperm_Y, freq_sperm_S, freq_sperm_YY, freq_sperm_SS, freq_sperm_YS])
					sex = egg + sperm
					if sex == 11:
						#XX female
						Num_XX += 1
						num_female += 1
						egg_X += RF_XX*6
					elif sex == 12 or sex == 21:
						#XXY female
						Num_XXY += 1
						num_female += 1
						egg_X += RF_XXY*3
						egg_XY += RF_XXY*3
					elif sex == 13 or sex == 31:
						#XXS female
						Num_XXS += 1
						num_female += 1
						egg_X += RF_XXS*3
						egg_XS += RF_XXS*3
					elif sex == 22:
						#XXYY female
						Num_XXYY += 1
						num_female += 1
						egg_XY += RF_XXYY*6
					elif sex == 33:
						#XXSS female
						Num_XXSS += 1
						num_female += 1
						egg_XS += RF_XXSS*6
					elif sex == 23 or sex == 32:
						#XXYS female
						Num_XXYS += 1
						num_female += 1
						egg_XY += RF_XXYS*3
						egg_XS += RF_XXYS*3
					elif sex == 14:
						#XY male
						Num_XY += 1
						num_male += 1
						sperm_X += RF_XY*3
						sperm_Y += RF_XY*3
					elif sex == 15:
						#XS male
						Num_XS += 1
						num_male += 1
						sperm_X += RF_XS*3
						sperm_S += RF_XS*3
					elif sex == 24 or sex == 16:
						#XYY male
						Num_XYY += 1
						num_male += 1
						sperm_XY += RF_XYY*2
						sperm_Y += RF_XYY*2
						sperm_X += RF_XYY*1
						sperm_YY += RF_XYY*1
					elif sex == 35 or sex == 17:
						#XSS male
						Num_XSS += 1
						num_male += 1
						sperm_XS += RF_XSS*2
						sperm_S += RF_XSS*2
						sperm_X += RF_XSS*1
						sperm_SS += RF_XSS*1
					elif sex == 25 or sex == 34 or sex == 18:
						#XYS male
						Num_XYS += 1
						num_male += 1
						sperm_XY += RF_XYS*1
						sperm_S += RF_XYS*1
						sperm_XS += RF_XYS*1
						sperm_Y += RF_XYS*1
						sperm_X += RF_XYS*1
						sperm_YS += RF_XYS*1
					elif sex == 26:
						#XYYY male
						Num_XYYY += 1
						num_male += 1
						sperm_XY += RF_XYYY*3
						sperm_YY += RF_XYYY*3
					elif sex == 28 or sex ==36:
						#XYYS male
						Num_XYYS += 1
						num_male += 1
						sperm_XY += RF_XYYS*2
						sperm_YS += RF_XYYS*2
						sperm_XS += RF_XYYS*1
						sperm_YY += RF_XYYS*1
					elif sex == 27 or sex ==38:
						#XYSS male
						Num_XYSS += 1
						num_male += 1
						sperm_XY += RF_XYSS*1
						sperm_SS += RF_XYSS*1
						sperm_XS += RF_XYSS*2
						sperm_YS += RF_XYSS*2
					else:
						#XSSS male 37
						Num_XSSS += 1
						num_male += 1
						sperm_XS += RF_XSSS*3
						sperm_SS += RF_XSSS*3
			
				Total_number_of_chromosomes_in_female = (2*Num_XX)+(3*Num_XXY)+(3*Num_XXS)+(4*Num_XXYY)+(4*Num_XXSS)+(4*Num_XXYS)
				Total_number_of_chromosomes_in_male = (2*Num_XY)+(2*Num_XS)+(3*Num_XYY)+(3*Num_XSS)+(3*Num_XYS)+(4*Num_XYYY)+(4*Num_XYYS)+(4*Num_XYSS)+(4*Num_XSSS)
				Total_number_of_X_in_female = (2*Num_XX)+(2*Num_XXY)+(2*Num_XXS)+(2*Num_XXYY)+(2*Num_XXSS)+(2*Num_XXYS)
				Total_number_of_Y_in_female = Num_XXY+(2*Num_XXYY)+Num_XXYS
				Total_number_of_S_in_female = Num_XXS+(2*Num_XXSS)+Num_XXYS
				Total_number_of_X_in_male = Num_XY+Num_XS+Num_XYY+Num_XSS+Num_XYS+Num_XYYY+Num_XYYS+Num_XYSS+Num_XSSS
				Total_number_of_Y_in_male = Num_XY+(2*Num_XYY)+Num_XYS+(3*Num_XYYY)+(2*Num_XYYS)+Num_XYSS
				Total_number_of_S_in_male = Num_XS+(2*Num_XSS)+Num_XYS+Num_XYYS+(2*Num_XYSS)+(3*Num_XSSS)
				
				Xf_freq = Total_number_of_X_in_female/Total_number_of_chromosomes_in_female
				Yf_freq = Total_number_of_Y_in_female/Total_number_of_chromosomes_in_female
				Sf_freq = Total_number_of_S_in_female/Total_number_of_chromosomes_in_female
				Xm_freq = Total_number_of_X_in_male/Total_number_of_chromosomes_in_male
				Ym_freq = Total_number_of_Y_in_male/Total_number_of_chromosomes_in_male
				Sm_freq = Total_number_of_S_in_male/Total_number_of_chromosomes_in_male
				S_freq = (Total_number_of_S_in_female+Total_number_of_S_in_male)/(Total_number_of_S_in_female+Total_number_of_S_in_male+Total_number_of_Y_in_female+Total_number_of_Y_in_male)

				
				Total_number_of_female_gametes = egg_X+egg_XY+egg_XS
				freq_egg_X = egg_X/Total_number_of_female_gametes
				freq_egg_XY = egg_XY/Total_number_of_female_gametes
				freq_egg_XS = egg_XS/Total_number_of_female_gametes		
				Total_number_of_male_gametes = sperm_X+sperm_XY+sperm_XS+sperm_Y+sperm_S+sperm_YY+sperm_SS+sperm_YS
				freq_sperm_X = sperm_X/Total_number_of_male_gametes
				freq_sperm_XY = sperm_XY/Total_number_of_male_gametes
				freq_sperm_XS = sperm_XS/Total_number_of_male_gametes
				freq_sperm_Y = sperm_Y/Total_number_of_male_gametes
				freq_sperm_S = sperm_S/Total_number_of_male_gametes
				freq_sperm_YY = sperm_YY/Total_number_of_male_gametes
				freq_sperm_SS = sperm_SS/Total_number_of_male_gametes
				freq_sperm_YS = sperm_YS/Total_number_of_male_gametes
				
				Num_1YSinF = Num_XXY+Num_XXS
				Num_2YSinF = Num_XXYY+Num_XXSS+Num_XXYS
				Num_1YSinM = Num_XY+Num_XS
				Num_2YSinM = Num_XYY+Num_XSS+Num_XYS
				Num_3YSinM = Num_XYYY+Num_XYYS+Num_XYSS+Num_XSSS
				
				N_1YSinF.append(Num_1YSinF)
				N_2YSinF.append(Num_2YSinF)
				N_1YSinM.append(Num_1YSinM)
				N_2YSinM.append(Num_2YSinM)
				N_3YSinM.append(Num_3YSinM)
				
				N_female.append(num_female)
				N_male.append(num_male)
				Xf_f.append(Xf_freq)
				Yf_f.append(Yf_freq)
				Sf_f.append(Sf_freq)
				Xm_f.append(Xm_freq)
				Ym_f.append(Ym_freq)
				Sm_f.append(Sm_freq)
				S_f.append(S_freq)
				if num_female == 0 or num_male == 0:
					break
			
			if S_freq == 0.0000:
				Sext += 1
			elif Ym_freq == 0.0000:
				Yext += 1
			else:
				YSco += 1
			MeanS = np.mean(S_f)
			Sfreq.append(MeanS)
		
		F_Sext = Sext/(Sext+Yext+YSco)
		F_Yext = Yext/(Sext+Yext+YSco)
		F_YSco = YSco/(Sext+Yext+YSco)
		F_Sey = np.mean(Sfreq)
		FitDif = Y-S
		Sex[i,j] = FitDif
		Yex[i,j] = F_Sey
		SYc[i,j] = F_YSco
		with open("result.txt", "a") as f:
			f.write("Y=\t")
			f.write(str(Y))
			f.write("\tS=\t")
			f.write(str(S))
			f.write("\tY-S=\t")
			f.write(str(FitDif))
			f.write("\tMeanSfreq=\t")
			f.write(str(F_Sey))
			f.write("\tSCoexistFreq=\t")
			f.write(str(F_YSco))
			f.write("\n")

#fig1, axA = plt.subplots()
#imA = axA.imshow(Sex, origin='lower', cmap="YlGnBu", vmin=-1.5, vmax=1.5)
#fig1.colorbar(imA)
#fig1.savefig('Y-S.png')

fig2, axB = plt.subplots()
imB = axB.imshow(Yex, origin='lower', cmap="YlGnBu", vmin=0, vmax=0.60)
plt.axis('off')
fig2.colorbar(imB)
fig2.savefig('S_freq.png')

#fig3, axC = plt.subplots()
#imC = axC.imshow(SYc, origin='lower', cmap="YlGnBu", vmin=0, vmax=1)
#plt.axis('off')
#fig3.colorbar(imC)
#fig3.savefig('Prob_YS_coext.png')
