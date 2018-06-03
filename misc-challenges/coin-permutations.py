# Author: Kirby Powell


#### GIVENS/GLOBALS ###

values = {"P": 1,
          "N": 5,
          "D": 10}

tester = "PPDDNNPDDPNN"
build_test = ""
ps = ["P", "P", "P", "P"]
ds = ["D", "D", "D", "D"]
ns = ["N", "N", "N", "N"]

clock = {0: "",
         1: "",
         2: "",
         3: "",
         4: "",
         5: "",
         6: "",
         7: "",
         8: "",
         9: "",
         10: "",
         11: ""}

### ACTUAL CODE BELOW ###

def validator(coins):
    start = 0

    while coins:
        current, coins = coins[0], coins[1:]

        if len(clock[(start + values[current]) % 12]) <= 0:
            clock[(start + values[current]) % 12] = current
            start = (start + values[current]) % 12
        else:
            return False

    return True


def permutations(builder, p, n, d, results):
    build1, build2, build3 = builder, builder, builder

    if len(builder) < 12:
        if len(p) > 0:
            build1 = (builder + "P")
            if validator(build1):
                permutations(build1, p[1:], n, d, results)
        if len(n) > 0:
            build2 = (builder + "N")
            if validator(build2):
                permutations(build2, p, n[1:], d, results)
        if len(d) > 0:
            build3 = (builder + "D")
            if validator(build3):
                permutations(build3, p, n, d[1:], results)
    else:
        if validator(builder):
            results.add(builder)

    return builder


if __name__ == "__main__":
    result = validator(tester)

    permResults = set()

    permutations(build_test, ps, ns, ds, permResults)
    print(permResults)

