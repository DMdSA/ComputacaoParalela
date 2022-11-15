# Parallel Computing - MEI - 22/23

## ToDo list

- Dá para meter flag de controlo na atualização dos pontos, em vez de nos clusters?
- Confirmar a diferença entre #I da versão sequencial e paralela (sequencial == sem primitivas OMP)
- Estudar CPI, Texec, Cache Misses, etc. => conclusões
- Confirmar: uso de simd, taskloop, taskwait
- static ({value}) -> estudar e ver se vale a pena (achoq vale)
- Só paralelizar partes de hot code?
- "ou paralelizamos nas folhas, ou na raiz"
- Granularidade e sincronizaçao paralelismo (pontos 3,4 dos slides de hoje)
- usar o time counter do openMP para comparar resultados internos (o do perf conta com mais coisas q não o código paralelo)
- 
- 

---

<br>

---

## Relatório

- Strong scalability(?) graph, reta +/- linear da comparação (Tseq / Tpar)
- Usar os resultados mais importantes nos gráficos/tabelas
- Na apresentação de valores, acordar com apenas 2/3 casas decimais e respeitar 0's
- 1 tabela/gráfico única c/ resultados mais importantes (costuma ser mais fácil de ler, avaliar e tirar conclusões)
- qual devia ser o speed up maximo que deveriamos obter com o paralelismo aplicado?