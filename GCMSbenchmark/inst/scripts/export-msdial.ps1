param(
    [Parameter(Mandatory=$true,Position=0)][string]$Install,
    [Parameter(Mandatory=$true,Position=1)][string]$Output,
    [Parameter(Mandatory=$true,Position=2,ValueFromRemainingArguments=$true)][string[]]$Files
)
$ErrorActionPreference = 'Stop'
$culture = [System.Globalization.CultureInfo]::InvariantCulture
[System.Threading.Thread]::CurrentThread.CurrentCulture = $culture
[System.Reflection.Assembly]::LoadFrom((Join-Path $Install 'Common.dll')) | Out-Null
[System.Reflection.Assembly]::LoadFrom((Join-Path $Install 'MsdialCore.dll')) | Out-Null
$utf8 = New-Object System.Text.UTF8Encoding $false
$writer = New-Object System.IO.StreamWriter $Output,$false,$utf8
try {
    $writer.WriteLine("sample`tid`trt`tarea`theight`tspectrum`tprofile")
    foreach ($file in $Files) {
        $version = 0
        $offsets = New-Object 'System.Collections.Generic.List[long]'
        $data = [CompMs.MsdialCore.Parser.MsdecResultsReader]::ReadMSDecResults(
            $file, [ref]$version, [ref]$offsets)
        if ($null -eq $data) { throw "Cannot read MS-DIAL deconvolution file $file" }
        $sample = [IO.Path]::GetFileNameWithoutExtension($file) -replace '_[0-9]+$',''
        foreach ($entry in $data) {
            $spectrum = ($entry.Spectrum | ForEach-Object {
                "$($_.Mass),$($_.Intensity)"
            }) -join ' '
            $profile = ($entry.ModelPeakChromatogram | ForEach-Object {
                "$($_.ChromXs.RT.Value),$($_.Intensity)"
            }) -join ' '
            $rt = $entry.ChromXs.RT.Value
            $writer.WriteLine(("$sample",$entry.ScanID,$rt,$entry.ModelPeakArea,
                $entry.ModelPeakHeight,$spectrum,$profile) -join "`t")
        }
    }
}
finally {
    $writer.Dispose()
}
