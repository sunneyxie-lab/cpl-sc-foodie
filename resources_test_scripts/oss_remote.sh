#!/bin/bash
# Sourced by upload_resources_test.sh and fetch_resources_test.sh. Configures rclone's
# "oss" remote from the OSS_* variables only, without writing a config file, and
# sets $remote to the fixtures' prefix.

for var in OSS_ENDPOINT OSS_BUCKET OSS_ACCESS_KEY_ID OSS_ACCESS_KEY_SECRET; do
  if [[ -z "${!var:-}" ]]; then
    echo "error: $var is not set" >&2
    exit 1
  fi
done

if ! command -v rclone > /dev/null; then
  echo "error: rclone not found on PATH. Install its static binary from https://downloads.rclone.org into ~/bin" >&2
  exit 1
fi

export RCLONE_CONFIG_OSS_TYPE=s3
# Generic S3 with path-style URLs: the lab's private OSS has no DNS for <bucket>.<endpoint>,
# and rclone's Alibaba provider always uses that form.
export RCLONE_CONFIG_OSS_PROVIDER=Other
export RCLONE_CONFIG_OSS_FORCE_PATH_STYLE=true
export RCLONE_CONFIG_OSS_ENDPOINT="$OSS_ENDPOINT"
export RCLONE_CONFIG_OSS_ACCESS_KEY_ID="$OSS_ACCESS_KEY_ID"
export RCLONE_CONFIG_OSS_SECRET_ACCESS_KEY="$OSS_ACCESS_KEY_SECRET"
# Without this rclone checks or creates the bucket first, which a key limited to the prefix can't do.
export RCLONE_CONFIG_OSS_NO_CHECK_BUCKET=true

# The lab cluster's OSS endpoint is signed by a private CA; OSS_CA_CERT points to its PEM file.
if [[ -n "${OSS_CA_CERT:-}" ]]; then
  if [[ ! -r "$OSS_CA_CERT" ]]; then
    echo "error: cannot read OSS_CA_CERT file $OSS_CA_CERT" >&2
    exit 1
  fi
  export RCLONE_CA_CERT="$OSS_CA_CERT"
fi

# The endpoint is on the lab's private network; a local proxy can't reach it.
export NO_PROXY="${NO_PROXY:+$NO_PROXY,}$OSS_ENDPOINT" no_proxy="${no_proxy:+$no_proxy,}$OSS_ENDPOINT"

# Fail within minutes instead of retrying for a long time when OSS is unreachable.
export RCLONE_CONTIMEOUT=15s RCLONE_TIMEOUT=60s RCLONE_RETRIES=2 RCLONE_LOW_LEVEL_RETRIES=3

remote="oss:$OSS_BUCKET/cpl-sc-foodie/resources_test/consensus_dedup"
