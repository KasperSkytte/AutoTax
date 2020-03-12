#!/usr/bin/env bats
load /autotax/autotax_functions.sh

@test "variable set: VERSION" {
  [ ${VERSION} ]
}
@test "variable set: silva_db" {
  [ ${silva_db} ]
}
@test "variable set: silva_udb" {
  [ ${silva_udb} ]
}
@test "variable set: typestrains_db" {
  [ ${typestrains_db} ]
}
@test "variable set: typestrains_udb" {
  [ ${typestrains_udb} ]
}
@test "variable set: denovo_prefix" {
  [ ${denovo_prefix} ]
}
@test "variable set: MAX_THREADS" {
  [ ${MAX_THREADS} ]
}
