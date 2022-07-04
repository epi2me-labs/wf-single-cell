/* Demonstration of how to use ArgumentParser and some "tests"
 *
 * Note: this needs to be moved to the top-level to actually run
 *       under nextflow run due to how imports magically work
 *
 */

nextflow.enable.dsl = 2

def arg_test(Map stuff) {
    def parser = new ArgumentParser(
        args:["arg1", "arg2"],
        kwargs:["kwarg1":1, "kwarg2":2],
        name:"arg_test")
    return parser.parse_args(stuff)
}


workflow {
    main:
        just_positionals = arg_test(
            ["arg1":1, "arg2":2])
        println(just_positionals)

        all_args = arg_test(
            ["arg1":1, "arg2":2, "kwarg1":4, "kwarg2":8])
        println(all_args)

        try {
        	fail_missing = arg_test(
            	["arg1":1, "kwarg1":4, "kwarg2":8])
            println(fail_missing)
        } catch (Exception ex) {
            println("Caught: ${ex.getCause()}")
        }

        try {
            fail_extra = arg_test(
                ["arg1":1, "arg2":2, "kwarg3":12, "kwarg2":8])
            println(fail_extra)
        } catch (Exception ex) {
            println("Caught: ${ex.getCause()}")
        }
    
}
