# -*- coding: utf-8 -*-
"""calculo_de_Eto.py
ENG 340 - Irrigação
"""

print('o objetivo desse script é automatizar o cálculo da equação de Pennman-Monteith FAO')

print('antes de efetivamente chegar na equação é preciso fazer uma série de passos para que cheguemos ao resultado final')
print('lembrando que muitos dados tem que ser retirados da estação metereológica, pois algumas variáveis não são passíveis de serem calculadas aqui')
print('ATENÇÃO: SEMPRE USE PONTO EM VEZ DE VÍRGULA PARA separar os números se não vai travar tudo e não vai dar certo')
print('Feito por Luciano M.D. Filho para a matéria de ENG 340 2024/1')

print('Cálculo de delta “∆” é a declividade da curva de pressão de vapor em relação à temperatura (Kpa°C^-1)')
# T é a temperatura media do ar em °C
import numpy as np
import math as math
J = float(input(print('Insira o valor de J (dias Julianos):   ')))
T = float(input(print('Insira o valor de T:')))
def calcular_delta(T):
  exp = (17.27*18.6)/(18.6 + 237.3)
  numerador = 4098*(0.6108*2.71828**exp)
  denominador = (18.6 + 237.3)**2
  delta = numerador / denominador
  return delta
delta1 = calcular_delta(T)
print(f"O valor de delta para T = {T} é {delta1} kPa°C^-1")

# calculando a Pressão atmosférica (Patm) em função de z (altitude do local)
z = float(input(print('Insira o valor de z (altitude):   ')))
def calcular_patm(z):
  ex = ((293-0.0065*z)/(293))**5.26
  funcao = 101.3 * ex
  Patm = funcao
  return Patm
patm = calcular_patm(z)
print(f"O valor de Patm para z = {z} é {patm} kPa")
# achando Y (constante psicometrica) essa é bem de boa depende do valor da Patm
y = 0.000665 * patm
print(f"O valor de y é {y} kPa°C^-1")

# calculo do deficit de saturação que é a diferença de es (pressão de vapor) e ea (pressão atual do vapor), claro antes iremos calcular es e ea separadamente
t = float(input('Insira o valor da Temperatura(T): '))
def calcular_es(t):
  exp = 2.71828**((17.27*t)/(t + 237.3))
  func = 0.6108 * exp
  es = func
  return es
es = calcular_es(t)
print(f"O valor de es é {es} kPa")
# calculo de ea
UR = float(input('Insira o valor da UR (em decimal, pois a UR normalmente é dada porcentagem e usar ponto em vez de vírgula):   '))
ea = (es * UR)
print(f"O valor de ea é {ea} kPa")

print('Estimativa do saldo de radiação (Rn), sendo que Rn = Rns - Rnl')
print('Para essa parte é importante tomar cuidado com as variáveis, vai ser pedido bastante coisa')
# para o calculo de Rn deve ser estimado antes Rns (saldo de radiação de ondas
# curtas (MJm-2dia-1)) e “Rnl” é o saldo de radiação  de ondas longas (MJm-2dia-1).
lat = float(input(print('Insira o valor da Latitude( valor em radiano):   ')))
Rs = float(input(print('Insira o valor de Rs (valores dadas pela estação) (em MJm-2dia-1):   ')))
def calcular_Rs(Rs):
  Rns = 0.77*Rs
  funcao = Rns
  return funcao
Rns = calcular_Rs(Rs)
print(f"O valor de Rns é {Rns} MJm-2dia-1")
print('Calculo de Rnl')
rnl_max = float(input(print('Insira o valor da Temperatura máxima do dia:  ')))
rnl_min = float(input(print('Insira o valor da Temperatura mínima do dia:  ')))
dr = 1 + 0.033*math.cos((2*math.pi/365)*J)
print(f"O valor de dr é {dr}")
declinacao_solar = 0.409*math.sin((2*math.pi/365*J-1.39))
print(f"O valor de declinacao_solar é {declinacao_solar}")
xx = (1-(math.tan(lat)**2)*(math.tan(declinacao_solar)**2))
if xx <= 0:
  xx == 0.00001
else:
  xx = xx
print(f"O valor de X é {xx}")
angulo_horario = math.pi/2 - math.atan((-math.tan(lat)*math.tan(declinacao_solar))/(xx**0.5))
Ra = (118.08/math.pi)*dr*(angulo_horario*math.sin(lat)*math.sin(declinacao_solar)+math.cos(lat)*math.cos(declinacao_solar)*math.sin(angulo_horario))
Rso = (0.75 + 2*(10**-5)*z)*Ra
print(f"O valor de Ra é {Ra} MJm-2dia-1")
print(f"O valor de Rso é {Rso} MJm-2dia-1")

def calcular_rnl(rnl_max, rnl_min, ea, Rs, Rso):
  parte0 = 0.000000004903
  parte1 = ((rnl_max+273.16)**4 + (rnl_min+273.16)**4)/2
  parte2 = (0.34-0.14*math.sqrt(ea))
  parte3 = 1.35*(Rs/Rso) - 0.35
  Rnl = parte0*parte1*parte2*parte3
  return Rnl
Rnl = calcular_rnl(rnl_max, rnl_min, ea, Rs, Rso)
Rn = Rns - Rnl
print(f"O valor de Rn é {Rn} MJm^-2*dia^-1")

print('Finalmente o calculo da Evapo de referência pelo método de Penman-Monteith')
delta_final = delta1
patm_final = patm
y_final = y
es_final = es
ea_final = ea
# na estimativa de calor no solo, a maioria das estações não vai ter o valor de G pode-se considerá-lo como sendo igual a 0, quando não houver medições disponíveis (ALLEN et al.,1998).
G = float(input(print('insira o valor de G(fluxo total diário de calor no solo): ')))
u2 = float(input(print('insira o valor de U2(velocidade do vento a 2 m de altura):  ')))
temperatura_final = float(input(print('Insira "novamente" a Temperatura:   ')))
Eto_numerador = 0.408*delta_final*(Rn - G) + (y_final*900*u2*(es_final - ea_final))/(temperatura_final+273)
Eto_denominador =  delta_final + y_final*(1+0.34*u2)
Eto_final = Eto_numerador/Eto_denominador
resultado = round(Eto_final, 3)
print(f"O valor de Eto, pelo método de Penman-Monteith, é  {resultado} mm/dia")
while True:
  exit = input('Digite OK para finalizar o app: ')
  if exit.upper() == 'OK':
    print('Obrigado por Utilizar')
    break
  else:
    print('Entrada inválida. Tente novamente')