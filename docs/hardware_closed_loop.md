# Integração de hardware — Closed-Loop (Arduino / NI-DAQ)

> ⚠️ **Antes de conectar qualquer hardware real**: leia
> `preditor/README_preditor.md`. Em 13/09/2026 o modelo de estado
> comportamental usado no braço de gating tem AUC-ROC 0.618 — não é preciso
> o bastante pra controlar um laser de optogenética com segurança. A aba
> "Closed-Loop / Replay" do programa (seção Análise de Dados) é um **replay offline**: mostra o que
> dispararia, não aciona nada de verdade. Este documento descreve ONDE
> plugar o hardware quando o modelo estiver validado, não recomenda ligá-lo
> hoje.

## Onde fica o ponto de integração

Um único ponto no código decide quando disparar: `dispara_ttl(t_seg)` em
[`preditor/prever_pac_tempo_real.py`](../preditor/prever_pac_tempo_real.py).
Hoje o corpo dessa função só imprime no console
(`>>> [TTL] Disparo de pulso em t=...`). Pra ligar hardware de verdade,
troque o CORPO dessa função — nada mais no pipeline precisa mudar, já que
`modulo_replay()` e a aba de replay do programa já chamam `dispara_ttl()` no momento
certo (`P(PAC) >= limiar`).

## Opção 1 — NI-DAQ (National Instruments, via `nidaqmx`)

Requer a biblioteca `nidaqmx` (`pip install nidaqmx`) e o driver NI-DAQmx
instalado (Windows). Canal digital de saída (ex.: porta 0, linha 0) acionado
em nível alto por ~10ms:

```python
import nidaqmx
import time

def dispara_ttl(t_seg):
    with nidaqmx.Task() as task:
        task.do_channels.add_do_chan("Dev1/port0/line0")
        task.write(True)
        time.sleep(0.01)   # 10 ms -- ajustar conforme o laser/driver
        task.write(False)
```

- `"Dev1/port0/line0"`: identificador do dispositivo/porta/linha, visível no
  NI MAX (Measurement & Automation Explorer). Ajustar pro nome real do seu
  DAQ.
- Duração do pulso (10ms aqui) deve casar com a especificação do driver do
  laser/LED de optogenética — confirmar com o datasheet do equipamento antes
  de qualquer teste com animal.

## Opção 2 — Arduino via porta serial

Requer a biblioteca `pyserial` (`pip install pyserial`). O Arduino escuta a
porta serial e aciona um pino digital ao receber `b'H'`/`b'L'` (exemplo de
sketch mínimo abaixo).

```python
import serial
import time

_ser = None

def dispara_ttl(t_seg):
    global _ser
    if _ser is None:
        _ser = serial.Serial('COM3', 9600)  # ajustar porta COM
    _ser.write(b'H')
    time.sleep(0.01)
    _ser.write(b'L')
```

Sketch Arduino correspondente (pino 13 como saída TTL):

```cpp
const int PINO_TTL = 13;

void setup() {
  pinMode(PINO_TTL, OUTPUT);
  Serial.begin(9600);
}

void loop() {
  if (Serial.available() > 0) {
    char c = Serial.read();
    if (c == 'H') digitalWrite(PINO_TTL, HIGH);
    if (c == 'L') digitalWrite(PINO_TTL, LOW);
  }
}
```

- Ajustar `'COM3'` pra porta real (ver Gerenciador de Dispositivos no
  Windows).
- Abrir a conexão serial UMA vez (fora do loop de replay) evita o atraso de
  reabrir a porta a cada disparo — por isso `_ser` é global/lazy acima.

## Checklist antes de ligar em um animal de verdade

1. Modelo revalidado com AUC-ROC e precisão/recall adequados ao protocolo
   (ver limitação documentada em `preditor/README_preditor.md`).
2. Pulso testado em bancada (osciloscópio ou LED de teste) confirmando
   largura e nível de tensão antes de conectar ao driver do laser.
3. Latência ponta-a-ponta (janela de LFP → decisão do modelo → pulso TTL)
   medida e compatível com a janela biológica de interesse.
4. Botão de parada de emergência independente do software (corte físico do
   circuito TTL), nunca dependente só do programa/Python continuar rodando.
