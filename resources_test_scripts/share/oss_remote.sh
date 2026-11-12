#!/bin/bash
# Sourced after the operation-specific credential file. Keep configuration in the
# environment only; never read or write an rclone config containing credentials.

for var in OSS_ENDPOINT OSS_BUCKET OSS_ACCESS_KEY_ID OSS_ACCESS_KEY_SECRET; do
  if [[ -z "${!var:-}" ]]; then
    printf 'error: %s is not set in the required credential file\n' "$var" >&2
    exit 1
  fi
done

if ! command -v rclone > /dev/null; then
  printf '%s\n' 'error: rclone not found on PATH. Install its static binary from https://downloads.rclone.org into ~/bin' >&2
  exit 1
fi

export RCLONE_CONFIG=/dev/null
export RCLONE_CONFIG_OSS_TYPE=s3
# The private endpoint has no DNS for <bucket>.<endpoint>; use generic path-style S3.
export RCLONE_CONFIG_OSS_PROVIDER=Other
export RCLONE_CONFIG_OSS_FORCE_PATH_STYLE=true
export RCLONE_CONFIG_OSS_ENDPOINT="$OSS_ENDPOINT"
export RCLONE_CONFIG_OSS_ACCESS_KEY_ID="$OSS_ACCESS_KEY_ID"
export RCLONE_CONFIG_OSS_SECRET_ACCESS_KEY="$OSS_ACCESS_KEY_SECRET"
# Prefix-limited keys cannot check or create the bucket.
export RCLONE_CONFIG_OSS_NO_CHECK_BUCKET=true

# Verify TLS even if the caller previously configured insecure rclone access.
export RCLONE_NO_CHECK_CERTIFICATE=false
unset RCLONE_CA_CERT
if [[ -n "${OSS_CA_CERT:-}" ]]; then
  if [[ ! -r "$OSS_CA_CERT" ]]; then
    printf 'error: cannot read OSS_CA_CERT file %s\n' "$OSS_CA_CERT" >&2
    exit 1
  fi
  export RCLONE_CA_CERT="$OSS_CA_CERT"
fi

# NO_PROXY takes a hostname, not an endpoint URL (which may include a port/path).
oss_endpoint_host=$(python3 - "$OSS_ENDPOINT" <<'PY'
import sys
from urllib.parse import urlsplit

endpoint = sys.argv[1]
host = urlsplit(endpoint if "://" in endpoint else "//" + endpoint).hostname
if not host:
    raise SystemExit("error: OSS_ENDPOINT must contain a hostname")
print(host)
PY
) || exit 1
export NO_PROXY="${NO_PROXY:+$NO_PROXY,}$oss_endpoint_host"
export no_proxy="${no_proxy:+$no_proxy,}$oss_endpoint_host"
unset oss_endpoint_host

# Bound connection failures instead of retrying indefinitely.
export RCLONE_CONTIMEOUT=15s RCLONE_TIMEOUT=60s RCLONE_RETRIES=2 RCLONE_LOW_LEVEL_RETRIES=3
