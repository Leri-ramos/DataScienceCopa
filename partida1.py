# -*- coding: utf-8 -*-
"""
Created on Tue Nov 22 10:52:28 2022

@author: user
"""

import streamlit as st  
import pandas as pd
import numpy as np
import random 
from scipy.stats import poisson 

st.set_page_config(
    page_title = 'Predições de Jogos da Copa do Mundo',
    page_icon = '⚽',
)

selecoes = pd.read_excel('DadosCopaDoMundoQatar2022.xlsx', sheet_name ='selecoes', index_col = 0)

fifa = selecoes['PontosRankingFIFA']
a, b = min(fifa), max(fifa) 
fa, fb = 0.15, 1 
b1 = (fb - fa)/(b-a) 
b0 = fb - b*b1
forca = b0 + b1*fifa 

def Resultado(gols1, gols2):
    if gols1 > gols2:
        res = 'V'
    if gols1 < gols2:
        res = 'D' 
    if gols1 == gols2:
        res = 'E'       
    return res

def MediasPoisson(selecao1, selecao2):
    forca1 = forca[selecao1]
    forca2 = forca[selecao2] 
    mgols = 2.75
    l1 = mgols*forca1/(forca1 + forca2)
    l2 = mgols*forca2/(forca1 + forca2)
    return [l1, l2]
    
def Distribuicao(media):
    probs = []
    for i in range(7):
        probs.append(poisson.pmf(i,media))
    probs.append(1-sum(probs))
    return pd.Series(probs, index = ['0', '1', '2', '3', '4', '5', '6', '7+'])

def ProbabilidadesPartida(selecao1, selecao2):
    l1, l2 = MediasPoisson(selecao1, selecao2)
    d1, d2 = Distribuicao(l1), Distribuicao(l2)  
    matriz = np.outer(d1, d2)    #   Monta a matriz de probabilidades

    vitoria = np.tril(matriz).sum()-np.trace(matriz)    #Soma a triangulo inferior
    derrota = np.triu(matriz).sum()-np.trace(matriz)    #Soma a triangulo superior
    probs = np.around([vitoria, 1-(vitoria+derrota), derrota], 3)
    probsp = [f'{100*i:.1f}%' for i in probs]

    nomes = ['0', '1', '2', '3', '4', '5', '6', '7+']
    matriz = pd.DataFrame(matriz, columns = nomes, index = nomes)
    matriz.index = pd.MultiIndex.from_product([[selecao1], matriz.index])
    matriz.columns = pd.MultiIndex.from_product([[selecao2], matriz.columns]) 
    output = {'seleção1': selecao1, 'seleção2': selecao2, 
             'f1': forca[selecao1], 'f2': forca[selecao2], 
             'media1': l1, 'media2': l2, 
             'probabilidades': probsp, 'matriz': matriz}
    return output

def Pontos(gols1, gols2):
    rst = Resultado(gols1, gols2)
    if rst == 'V':
        pontos1, pontos2 = 3, 0
    if rst == 'E':
        pontos1, pontos2 = 1, 1
    if rst == 'D':
        pontos1, pontos2 = 0, 3
    return pontos1, pontos2, rst


def Jogo(selecao1, selecao2):
    l1, l2 = MediasPoisson(selecao1, selecao2)
    gols1 = int(np.random.poisson(lam = l1, size = 1))
    gols2 = int(np.random.poisson(lam = l2, size = 1))
    saldo1 = gols1 - gols2
    saldo2 = -saldo1
    pontos1, pontos2, result = Pontos(gols1, gols2)
    placar = '{}x{}'.format(gols1, gols2)
    return [gols1, gols2, saldo1, saldo2, pontos1, pontos2, result, placar]

def JogosGrupo(dados, grupo): 
    
    times = list(dados.loc[dados['Grupo']==grupo].index)

    time1, time2, time3, time4 = times  # nome dos times 

    pt1, pt2, pt3, pt4 = 0, 0, 0, 0  # pontos
    gp1, gp2, gp3, gp4 = 0, 0, 0, 0  #gols pro
    sg1, sg2, sg3, sg4 = 0, 0, 0, 0 #saldo de gols

    
    
    jogo1 = Jogo(time1, time2)
    jogo2 = Jogo(time3, time4)

    jogo3 = Jogo(time1, time3)
    jogo4 = Jogo(time2, time4)

    jogo5 = Jogo(time1, time4)
    jogo6 = Jogo(time2, time3)

    gp1, gp2, sg1, sg2, pt1, pt2 = gp1 + jogo1[0], gp2 + jogo1[1], sg1 + jogo1[2], sg2 + jogo1[3], pt1 + jogo1[4], pt2 + jogo1[5]
    gp3, gp4, sg3, sg4, pt3, pt4 = gp3 + jogo2[0], gp4 + jogo2[1], sg3 + jogo2[2], sg4 + jogo2[3], pt3 + jogo2[4], pt4 + jogo2[5]
    gp1, gp3, sg1, sg3, pt1, pt3 = gp1 + jogo3[0], gp3 + jogo3[1], sg1 + jogo3[2], sg3 + jogo3[3], pt1 + jogo3[4], pt3 + jogo3[5]
    gp2, gp4, sg2, sg4, pt2, pt4 = gp2 + jogo4[0], gp4 + jogo4[1], sg2 + jogo4[2], sg4 + jogo4[3], pt2 + jogo4[4], pt4 + jogo4[5]
    gp1, gp4, sg1, sg4, pt1, pt4 = gp1 + jogo5[0], gp4 + jogo5[1], sg1 + jogo5[2], sg4 + jogo5[3], pt1 + jogo5[4], pt4 + jogo5[5]
    gp2, gp3, sg2, sg3, pt2, pt3 = gp2 + jogo6[0], gp3 + jogo6[1], sg2 + jogo6[2], sg3 + jogo6[3], pt2 + jogo6[4], pt3 + jogo6[5]

    partidas = [ time1 + ' x ' + time2, 
                 time3 + ' x ' + time4,
                 time1 + ' x ' + time3, 
                 time2 + ' x ' + time4,
                 time1 + ' x ' + time4,
                 time2 + ' x ' + time3 ]

    resultados = [ jogo1[6], jogo2[6], jogo3[6], jogo4[6], jogo5[6], jogo6[6] ]
    placares = [ jogo1[-1], jogo2[-1], jogo3[-1], jogo4[-1], jogo5[-1], jogo6[-1] ] 
    cols = ['Pontos', 'Saldo de Gols', 'Gols Pró']
    tab = pd.DataFrame([[pt1, pt2, pt3, pt4], [sg1, sg2, sg3, sg4], [gp1, gp2, gp3, gp4]], index = cols, columns = times).transpose()
    
    tab = tab.sort_values(['Pontos', 'Saldo de Gols', 'Gols Pró'], ascending = False)
    tab['Posição'] = [1, 2, 3, 4]

    jogos = pd.DataFrame([partidas, placares, resultados]).transpose()
    jogos.columns = ['Partida', 'Placar', 'Resultado']

    return [tab, jogos]

def JogoMataMata(selecao1, selecao2): 
    jogo = Jogo(selecao1, selecao2)
    resultado = jogo[6]
    if resultado == 'V':
        return selecao1
    elif resultado == 'D':
        return selecao2
    else:
        return random.sample([selecao1, selecao2], 1)[0]

# Segunda fase da competição 

def SimulaCopa(selecoes):
  #Guarda resultados
    cols = ['1st', '2nd', '3rd', '4th', 'Oitavas', 'Quartas', 'Semis', 'Final', 'Campeão']
    n = selecoes.shape[0]
    m = len(cols)
    aux = np.array(np.zeros(n*m).reshape(n, m))
    info = pd.DataFrame(aux, columns = cols, index = selecoes.index) 
    info = info.astype(int)
    info
      #simular toda primeira fase
    top16 = []  
    for i in list('ABCDEFGH'):
        a = JogosGrupo(selecoes, grupo = i)[0] 
        top16 += a.index[:2].tolist()
        anomes = a.index.to_list() 
        info.at[anomes[0], '1st'] = 1
        info.at[anomes[1], '2nd'] = 1
        info.at[anomes[2], '3rd'] = 1
        info.at[anomes[3], '4th'] = 1


      #tabela da seleção - Simula as oitavas
    qf1 = JogoMataMata(top16[0], top16[3])   #A1 x B2
    qf2 = JogoMataMata(top16[2], top16[1])   #B1 x A2
    qf3 = JogoMataMata(top16[4], top16[7])   #C1 x D2 
    qf4 = JogoMataMata(top16[6], top16[5])   #D1 x C2
    qf5 = JogoMataMata(top16[8], top16[11])  #E1 x F2
    qf6 = JogoMataMata(top16[10], top16[9])  #F1 x E2
    qf7 = JogoMataMata(top16[12], top16[15]) #G1 x H2
    qf8 = JogoMataMata(top16[14], top16[13]) #H1 x G2

    top8 = [qf1, qf2, qf3, qf4, qf5, qf6, qf7, qf8]

    # simula as quartas
    sf1 = JogoMataMata(qf1, qf3)
    sf2 = JogoMataMata(qf2, qf4) 
    sf3 = JogoMataMata(qf5, qf7) 
    sf4 = JogoMataMata(qf6, qf8) 

    top4 = [sf1, sf2, sf3, sf4]
    #simula os semifinalistas 
    f1 = JogoMataMata(sf1, sf3) 
    f2 = JogoMataMata(sf2, sf4) 
    top2 = [f1, f2]

    #campeão / final 

    top1 = JogoMataMata(f1, f2)

    #guardar as informações

    info.at[top16, 'Oitavas'] = 1
    info.at[top8, 'Quartas'] = 1
    info.at[top4, 'Semis'] = 1
    info.at[top2, 'Final'] = 1
    info.at[top1, 'Campeão'] = 1

    return info



    

def SimulacaoTotal(dados, S = 1000): 
    print('IA: "Iniciando simulação..."')
    info = SimulaCopa(dados)
    for i in range(S-1):
        info += SimulaCopa(dados)
        if (i+2)%(S/10) == 0:
            print('IA: "Simulações de Copas do Mundo: {:.0f}% completas'.format(100*((i+2)/S)))    
    print('IA: "Simulação Finalizada!"')
    return info.sort_values(by = 'Campeão', ascending = False)/S

######## COMEÇO DO APP


st.markdown("# 🏆 Copa do Mundo Qatar 2022") 

st.markdown("## ⚽ Probabilidades das Partidas")
st.markdown('---')

listaselecoes1 = selecoes.index.tolist()  
listaselecoes1.sort()
listaselecoes2 = listaselecoes1.copy()
#botoes de escolha, se escolher selecao1, revome
j1, j2 = st.columns(2)
selecao1 = j1.selectbox('Escolha a primeira Seleção', listaselecoes1) 
listaselecoes2.remove(selecao1)
selecao2 = j2.selectbox('Escolha a segunda Seleção', listaselecoes2, index = 1)
st.markdown('---')

jogo = ProbabilidadesPartida(selecao1, selecao2)
prob = jogo['probabilidades']
matriz = jogo['matriz']

col1, col2, col3, col4, col5 = st.columns(5)
col1.image(selecoes.loc[selecao1, 'LinkBandeiraGrande'])  
col2.metric(selecao1, prob[0])
col3.metric('Empate', prob[1])
col4.metric(selecao2, prob[2]) 
col5.image(selecoes.loc[selecao2, 'LinkBandeiraGrande'])

st.markdown('---')
st.markdown("## 📊 Probabilidades dos Placares") 

def aux(x):
	return f'{str(round(100*x,1))}%'
st.table(matriz.applymap(aux)) #vai aplicar em todos os valores das matrizes


st.markdown('---')
st.markdown("## 🌍 Probabilidades dos Jogos da Copa") 

jogoscopa = pd.read_excel('dados/outputEstimativasJogosCopa.xlsx', index_col = 0)
st.table(jogoscopa[['grupo', 'seleção1', 'seleção2', 'Vitória', 'Empate', 'Derrota']])


st.markdown('---')
st.markdown('Divirta-se')

#bandeira1, nome1, prob, empate, prob, nome2, bandeira2
#matriz de probabilidades do jogo
#placar mais provável

