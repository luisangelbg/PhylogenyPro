# Imprime manual-completo.html a PDF manejando Chrome por el protocolo de DevTools.
# Se hace asi, y no con --print-to-pdf, porque ese imprime antes de que terminen de cargar las
# tipografias de Google y antes de que paginar.js reparta las hojas: el PDF salia en Palatino,
# Georgia y Segoe UI, que son los tipos de reserva.
param(
  [string]$Url  = 'http://localhost:9801/manual/es/manual-completo.html',
  [string]$Out  = 'C:\Users\luisa\AppData\Local\Temp\claude\C--Users-luisa-Documents\160abd69-02c8-4b15-92b7-c2720dda96cc\scratchpad\PhylogenyPro-manual.pdf',
  [string]$Prof = 'C:\Users\luisa\AppData\Local\Temp\claude\C--Users-luisa-Documents\160abd69-02c8-4b15-92b7-c2720dda96cc\scratchpad\chrome-cdp',
  [int]$Port    = 9333
)
$ErrorActionPreference = 'Stop'
$chrome = 'C:\Program Files\Google\Chrome\Application\chrome.exe'

$proc = Start-Process $chrome -PassThru -WindowStyle Hidden -ArgumentList @(
  '--headless=new', '--disable-gpu', '--no-first-run', '--no-default-browser-check',
  '--disable-extensions', '--hide-scrollbars', '--window-size=1000,1300',
  "--remote-debugging-port=$Port", "--user-data-dir=$Prof", $Url)

function Esperar-Endpoint {
  for ($i = 0; $i -lt 60; $i++) {
    try { return (Invoke-WebRequest "http://127.0.0.1:$Port/json/list" -UseBasicParsing -TimeoutSec 3).Content } catch { Start-Sleep -Milliseconds 500 }
  }
  throw 'Chrome no abrio el puerto de depuracion'
}
$lista = Esperar-Endpoint | ConvertFrom-Json
$tab = $lista | Where-Object { $_.type -eq 'page' } | Select-Object -First 1
if (-not $tab) { throw 'No hay pestana' }
"pestana: $($tab.url)"

$ws = New-Object System.Net.WebSockets.ClientWebSocket
$ct = [System.Threading.CancellationToken]::None
$ws.ConnectAsync([Uri]$tab.webSocketDebuggerUrl, $ct).Wait()

$script:id = 0
function Enviar($metodo, $params) {
  $script:id++
  $msg = @{ id = $script:id; method = $metodo }
  if ($params) { $msg.params = $params }
  $json = $msg | ConvertTo-Json -Depth 10 -Compress
  $bytes = [System.Text.Encoding]::UTF8.GetBytes($json)
  $seg = New-Object System.ArraySegment[byte] (,$bytes)
  $ws.SendAsync($seg, [System.Net.WebSockets.WebSocketMessageType]::Text, $true, $ct).Wait()
  # leer hasta encontrar la respuesta con nuestro id (los eventos se descartan)
  while ($true) {
    $sb = New-Object System.Text.StringBuilder
    do {
      $buf = New-Object byte[] 65536
      $s2 = New-Object System.ArraySegment[byte] (,$buf)
      $r = $ws.ReceiveAsync($s2, $ct); $r.Wait()
      [void]$sb.Append([System.Text.Encoding]::UTF8.GetString($buf, 0, $r.Result.Count))
    } while (-not $r.Result.EndOfMessage)
    $o = $sb.ToString() | ConvertFrom-Json
    if ($o.PSObject.Properties.Name -contains 'id' -and $o.id -eq $script:id) {
      if ($o.PSObject.Properties.Name -contains 'error') { throw "$metodo : $($o.error.message)" }
      return $o.result
    }
  }
}

[void](Enviar 'Page.enable' $null)
[void](Enviar 'Runtime.enable' $null)

# Esperar a que las tres tipografias esten cargadas Y a que el numero de hojas deje de cambiar.
$espera = @'
new Promise(async (listo) => {
  const t0 = Date.now();
  const fams = [['Cormorant','500 16px "Cormorant"'], ['Crimson Pro','400 16px "Crimson Pro"'], ['Jost','400 16px "Jost"']];
  const cuenta = () => document.querySelectorAll('.libro .hoja').length;
  let previo = -1, estable = 0;
  while (Date.now() - t0 < 180000) {
    await new Promise(r => setTimeout(r, 500));
    const n = cuenta();
    if (n > 0 && n === previo) estable++; else estable = 0;
    previo = n;
    const cargadas = fams.every(f => document.fonts.check(f[1]));
    if (n > 0 && estable >= 4 && cargadas) {
      return listo(JSON.stringify({ hojas: n, fuentes: fams.filter(f => document.fonts.check(f[1])).map(f => f[0]), estado: document.fonts.status }));
    }
  }
  listo(JSON.stringify({ hojas: cuenta(), fuentes: fams.filter(f => document.fonts.check(f[1])).map(f => f[0]), estado: document.fonts.status, agotado: true }));
})
'@
$r = Enviar 'Runtime.evaluate' @{ expression = $espera; awaitPromise = $true; returnByValue = $true; timeout = 200000 }
"estado: $($r.result.value)"

$pdf = Enviar 'Page.printToPDF' @{
  printBackground   = $true
  preferCSSPageSize = $true
  marginTop = 0; marginBottom = 0; marginLeft = 0; marginRight = 0
  paperWidth = 8.5; paperHeight = 11
  displayHeaderFooter = $false
  scale = 1
}
[System.IO.File]::WriteAllBytes($Out, [Convert]::FromBase64String($pdf.data))
"PDF: {0:N2} MB" -f ((Get-Item $Out).Length / 1MB)

$ws.CloseAsync([System.Net.WebSockets.WebSocketCloseStatus]::NormalClosure, 'fin', $ct).Wait()
try { $proc | Stop-Process -Force } catch { }
