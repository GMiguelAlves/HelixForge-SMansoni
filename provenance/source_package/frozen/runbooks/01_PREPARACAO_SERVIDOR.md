# Preparação do servidor Slurm

1. Copie o pacote e renomeie `env/server.env.template` para `env/server.env`.
2. Edite caminhos, partição, conta e QoS; carregue com `source env/server.env`.
3. Clone o HelixForge e faça checkout do commit registrado em `HF_HELIXFORGE_COMMIT`.
4. Instale Nextflow 25.10.7 no caminho `HF_NEXTFLOW_BIN` e carregue Java 21.
5. Garanta que nós de compute leiam pacote, dados, referência, resultados e workdir.
6. O perfil padrão é `apptainer,slurm`. Apptainer é experimental no contrato v1; execute primeiro um smoke/piloto no cluster. Para Conda, use `HF_PROFILES=conda,slurm` e valide os ambientes.

Não rode tarefas científicas pesadas no nó de login. O processo Nextflow permanece no nó de submissão e envia cada processo ao Slurm.
