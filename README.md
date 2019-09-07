A project that calculates a CPU heavy task, built to work faster through MPI implementation.
Both solutions work with Master-Slave methodology.

Static -   All tasks are split evenly between all processes, each process works on his batch.
           The master process is responsible for handling leftover tasks.
Dynamic -  The master process is responsible only for handling information between itself and the other processes.
           Each slave is given an initial task, when it completes the task, it requests another one from the master process until all                tasks are done.
