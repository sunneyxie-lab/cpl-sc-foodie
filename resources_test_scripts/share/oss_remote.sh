#!/bin/bash
# Sourced with upload or download; keep configuration in the environment only.
# Never read or write an rclone config containing credentials.

for var in "${!RCLONE_CONFIG_OSS_@}" "${!RCLONE_S3_@}"; do
  unset "$var"
done

if ! command -v rclone > /dev/null; then
  printf '%s\n' 'error: rclone not found on PATH. Install its static binary from https://downloads.rclone.org into ~/bin' >&2
  exit 1
fi

export RCLONE_CONFIG=/dev/null
export RCLONE_CONFIG_OSS_TYPE=s3
export RCLONE_CONFIG_OSS_PROVIDER=Alibaba
export RCLONE_CONFIG_OSS_FORCE_PATH_STYLE=false
export RCLONE_CONFIG_OSS_ENDPOINT=https://oss-cn-beijing.aliyuncs.com
export RCLONE_CONFIG_OSS_ENV_AUTH=false
if [[ "$1" == upload ]]; then
  for var in OSS_ACCESS_KEY_ID OSS_ACCESS_KEY_SECRET; do
    if [[ -z "${!var:-}" ]]; then
      printf 'error: %s is not set in the required credential file\n' "$var" >&2
      exit 1
    fi
  done
  export RCLONE_CONFIG_OSS_ACCESS_KEY_ID="$OSS_ACCESS_KEY_ID"
  export RCLONE_CONFIG_OSS_SECRET_ACCESS_KEY="$OSS_ACCESS_KEY_SECRET"
else
  export RCLONE_CONFIG_OSS_ACCESS_KEY_ID=
  export RCLONE_CONFIG_OSS_SECRET_ACCESS_KEY=
fi
# Prefix-limited keys cannot check or create the bucket.
export RCLONE_CONFIG_OSS_NO_CHECK_BUCKET=true

# Verify TLS even if the caller previously configured insecure rclone access.
export RCLONE_NO_CHECK_CERTIFICATE=false
unset RCLONE_CA_CERT

# Bypass proxies for both the endpoint and virtual-hosted bucket requests.
export NO_PROXY="${NO_PROXY:+$NO_PROXY,}oss-cn-beijing.aliyuncs.com,xielab-foodie.oss-cn-beijing.aliyuncs.com"
export no_proxy="${no_proxy:+$no_proxy,}oss-cn-beijing.aliyuncs.com,xielab-foodie.oss-cn-beijing.aliyuncs.com"

# Bound connection failures instead of retrying indefinitely.
export RCLONE_CONTIMEOUT=15s RCLONE_TIMEOUT=60s RCLONE_RETRIES=2 RCLONE_LOW_LEVEL_RETRIES=3
