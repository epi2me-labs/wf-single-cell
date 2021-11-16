
process ping {
    label "pysam"
    cpus 1
    input:
        val message
    script:
        hostname = InetAddress.getLocalHost().getHostName()
        opsys = System.properties['os.name'].toLowerCase()
        disable = params.disable_ping ? '--disable' : ''
    """
    ping.py \
        --hostname $hostname \
        --opsys "$opsys" \
        --session $workflow.sessionId \
        --message $message \
        $disable
    """
}