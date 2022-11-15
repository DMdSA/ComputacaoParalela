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




### context switching
What happens in context switching?
Context switching is the process of storing the state of a thread so that it can be restored to resume execution at a later point in time. Rapid context switching between threads is expensive in terms of CPU utilization. Each context switch takes the kernel about 5 μs (on average) to process.
<br>

Certain constructs have task scheduling points at defined locations within them. When a thread encounters a task scheduling point, 
it is allowed to suspend the current task and execute another (called task switching). It can then return to the original task and resume.

<br>

### page fault
A page fault is an interruption that occurs when a software program attempts to access a memory block not currently stored in the system's RAM. This exception tells the operating system to find the block in virtual memory so it can be sent from a device's storage (SSD or HD) to RAM

<br>

### granularity and parallel performance

One key to attaining good parallel performance is choosing the right granularity for the
application. Granularity is the amount of real work in the parallel task. If granularity is too fine,
then performance can suffer from communication overhead. If granularity is too coarse, then
performance can suffer from load imbalance. The goal is to determine the right granularity
(usually larger is better) for parallel tasks, while avoiding load imbalance and communication
overhead to achieve the best performance.

> https://www.intel.com/content/dam/develop/external/us/en/documents/1-3-appthr-granularity-and-parallel-performance-699424.pdf

<br>

The schedule (dynamic, 1) clause specifies that the scheduler distribute one
iteration (or chunk) at a time dynamically to each thread. Each worker thread processes one
iteration and then returns to the scheduler, and synchronizes to get another iteration. By
increasing the chunk size, we increase the work size for each task that is assigned to a thread
and therefore reduce the number of times each thread must synchronize with the scheduler.
While this approach can improve performance, one must bear in mind (as mentioned above) that
increasing the granularity too much can cause load imbalance.

<br>

A more appropriate work size (100) should be used to select the right granularity for the program.
Also, since the difference in the amount of work between consecutive tasks will be less severe
than the previous chunk size, a further elimination of parallel overhead can be accomplished by
using the static schedule rather than dynamic. The code below shows the change in the schedule
clause that will virtually the overhead from this code segment and produce the fastest overall
parallel performance.

<br>

