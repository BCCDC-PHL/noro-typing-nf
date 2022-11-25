
// possibly useful template  https://www.nextflow.io/docs/latest/process.html
// using storeDir will only run the process if the output files DO NOT exist under the specified path
// if they exist, process is skipped and output is passed as the existing files
process formatBlastDatabases {
  storeDir '/db/genomes'

  input:
  path species

  output:
  path "${dbName}.*"

  script:
  dbName = species.baseName
  """
  makeblastdb -dbtype nucl -in ${species} -out ${dbName}
  """
}